# A simple illustration to show the deficiences of PCA, SPCA, CCA, and their partial counterparts
# for causal dimension reduction.

library(pracma)
# ========================
# Supervised PCA & CCA Methods
# ========================

# Standard SPCA
supervised_pca <- function(X, Y) {
  cov_xy <- t(X) %*% Y
  omega <- cov_xy / sqrt(sum(cov_xy^2))
  return(omega)
}

# CCA with scalar Y
cca_scalar <- function(X, Y) {
  Sigma_xx <- cov(X)
  Sigma_xy <- cov(X, Y)[, 1]
  omega <- solve(Sigma_xx, Sigma_xy)
  omega <- omega / sqrt(sum(omega^2))
  return(omega)
}

# Partial SPCA (adjusting for confounders C)
partial_supervised_pca <- function(X, Y, C) {
  Sigma_xy <- cov(X, Y)[,1]
  Sigma_xc <- cov(X, C)
  Sigma_cc <- cov(C)
  Sigma_cy <- cov(C, Y)[,1]
  
  Sigma_xy_c <- Sigma_xy - Sigma_xc %*% solve(Sigma_cc, Sigma_cy)
  omega <- Sigma_xy_c / sqrt(sum(Sigma_xy_c^2))
  return(omega)
}

# Partial CCA (adjusting for confounders C)
partial_cca_scalar <- function(X, Y, C) {
  Sigma_xx <- cov(X)
  Sigma_xy <- cov(X, Y)[,1]
  Sigma_xc <- cov(X, C)
  Sigma_cc <- cov(C)
  Sigma_cx <- t(Sigma_xc)
  Sigma_cy <- cov(C, Y)[,1]
  
  Sigma_xx_c <- Sigma_xx - Sigma_xc %*% solve(Sigma_cc, Sigma_cx)
  Sigma_xy_c <- Sigma_xy - Sigma_xc %*% solve(Sigma_cc, Sigma_cy)
  
  omega <- solve(Sigma_xx_c, Sigma_xy_c)
  omega <- omega / sqrt(sum(omega^2))
  return(omega)
}

doPCA <- function(X, Y, C){
  
  # === Main function ===
  find_optimal_omega <- function(X, Y, C, n_grid = 100, n_iter = 200, lr = 1e-2, verbose = TRUE) {
    n <- nrow(X)
    p <- ncol(X)
    
    # Step 1: residualize Y ~ C to estimate E[Y | Z = z, C], then marginalize C
    Y_resid <- residuals(lm(Y ~ C))
    
    # Objective function: for a given omega, estimate Var(mu(z))
    estimate_mu_var <- function(omega_unit) {
      Z <- as.vector(X %*% omega_unit)
      
      # Nonparametric estimate of mu(z) = E[Y_resid | Z = z]
      model <- smooth.spline(Z, Y_resid)
      z_grid <- seq(min(Z), max(Z), length.out = n_grid)
      mu_z <- predict(model, x = z_grid)$y
      
      return(var(mu_z))
    }
    
    # Gradient approximation via finite difference
    estimate_gradient <- function(omega, eps = 1e-4) {
      grad <- numeric(length(omega))
      for (j in 1:length(omega)) {
        d <- rep(0, length(omega)); d[j] <- eps
        omega_plus <- normalize(omega + d)
        omega_minus <- normalize(omega - d)
        grad[j] <- (estimate_mu_var(omega_plus) - estimate_mu_var(omega_minus)) / (2 * eps)
      }
      return(grad)
    }
    
    # Step 2: Initialize and optimize
    omega <- normalize(rnorm(p))  # random unit vector
    
    for (iter in 1:n_iter) {
      grad <- estimate_gradient(omega)
      omega <- normalize(omega + lr * grad)
      if (verbose && iter %% 10 == 0) {
        cat(sprintf("Iter %d | Var(mu(z)) = %.5f\n", iter, estimate_mu_var(omega)))
      }
    }
    
    return(list(omega = omega, var_mu = estimate_mu_var(omega)))
  }
  
  # === Usage example ===
  # set.seed(42)
  # n <- 500; p <- 5; q <- 2
  # C <- matrix(rnorm(n * q), n, q)
  # X <- matrix(rnorm(n * p), n, p)
  # X[,1] <- 0.5 * C[,1] + rnorm(n)
  # Y <- X[,1] * sign(C[,1]) + rnorm(n)

  # result <- find_optimal_omega(X, Y, C)
  # result$omega  # optimal direction
  
}

# ========================
# Gram-Schmidt Orthogonalization
# ========================
orthogonalize <- function(vec, basis) {
  for (b in 1:ncol(basis)) {
    proj <- sum(vec * basis[, b]) * basis[, b]
    vec <- vec - proj
  }
  return(vec)
}

# ========================
# Multi-Component Extraction (Orthogonalized)
# ========================

extract_components <- function(X, Y, ncomp = 1, method = c("spca", "cca", "partial_spca", "partial_cca"), C = NULL) {
  if(ncomp > 1){
    stop("More than 1 component is currently not supported")
  }
  method <- match.arg(method)
  
  # Choose extraction function
  find_direction <- switch(method,
                           spca = supervised_pca,
                           cca = cca_scalar,
                           partial_spca = function(X, Y) partial_supervised_pca(X, Y, C),
                           partial_cca = function(X, Y) partial_cca_scalar(X, Y, C)
  )
  
  X <- scale(X)
  Y <- scale(Y)
  if (!is.null(C)) C <- scale(C)
  
  p <- ncol(X)
  W <- matrix(0, p, ncomp)   # loadings
  Z <- matrix(0, nrow(X), ncomp) # component scores
  
  for (k in 1:ncomp) {
    omega_k <- find_direction(X, Y)
    
    if (k > 1) {
      omega_k <- orthogonalize(omega_k, W[, 1:(k-1), drop = FALSE])
    }
    
    omega_k <- omega_k / sqrt(sum(omega_k^2))  # enforce unit norm
    W[, k] <- omega_k
    Z[, k] <- X %*% omega_k
  }
  
  colnames(W) <- paste0("Comp", 1:ncomp)
  colnames(Z) <- paste0("Comp", 1:ncomp)
  return(list(#scores = Z, 
              loadings = W)
         )
}
# # Y_threshold <- .5*X[,1] + 3*X[,2] + X[,1] * ifelse(C[,1] > 0, 5, -5) + rnorm(n)
# Y_threshold <- X[,1] + X[,2] + 2 * X[,1] * exp(C) + rnorm(n)
# ========================
# Example: causal setting with p = 3 and confounding
# ========================

set.seed(42)
n <- 200; p <- 3; q <- 1  # 10 exposures, 1 confounder

# Step 1: Generate confounder(s)
C <- matrix(rnorm(n * q), n, q)

# Step 2: Generate exposures X with some depending on C
X     <- matrix(rnorm(n * p), n, p)
X[,1] <- C[,1] + rnorm(n)     # Mild confounding + signal
X[,2] <- C[,1] + rnorm(n)    # Mild confounding + no signal
                                # X3 is just independent noise

# Step 3: Generate outcome Y as a function of X1 and C. 

Y_linear   <- ( X[,1] + C[,1] + rnorm(n) )
Y_threshold <- X[,1] + 10*X[,2]*C[,1] + rnorm(n)

Yall <- cbind(Y_linear, Y_threshold)

loadings_result <- vector(mode= "list", 2L)

for(i in 1:2){
  Y <- Yall[, i]
  # Step 4: Extract components
  pca <- prcomp(X, center = TRUE)$rotation[,1]
  result_spca <- extract_components(X, Y, ncomp = 1, method = "spca")
  result_cca <- extract_components(X, Y, ncomp = 1, method = "cca")
  result_partial_spca <- extract_components(X, Y, ncomp = 1, method = "partial_spca", C = C)
  result_partial_cca <- extract_components(X, Y, ncomp = 1, method = "partial_cca", C = C)
  
  #  Create results matrix
  L <- cbind(pca,
            result_spca$loadings,
            result_cca$loadings,
            result_partial_spca$loadings,
            result_partial_cca$loadings)
  
  # Example where confounding is specified correctly, 
  # PCA, SPCA, and CCA fail to downweight X2, but P-SPCA and P-CCA do. 
  colnames(L) <- c("PCA", "SPCA", "CCA", "P-SPCA", "P-CCA")
  
  loadings_result[[i]] <- L
}

names(loadings_result) <- c("Linear confounding", "Threshold/nonlinear confounding
                            ")

# $`Linear confounding`
# PCA        SPCA          CCA      P-SPCA         P-CCA
# [1,] -0.46298970  0.78352545  0.926794400  0.99416291  0.9964176043
# [2,] -0.88602129  0.62080443  0.375515356 -0.08918402 -0.0845682816
# [3,]  0.02463362 -0.02626275 -0.006352747  0.06071508 -0.0004043893
# 
# $`Threshold/nonlinear confounding\n                            `
# PCA      SPCA       CCA     P-SPCA     P-CCA
# [1,] -0.46298970 0.5212327 0.4795493 0.33160639 0.4994578
# [2,] -0.88602129 0.3531516 0.1675463 0.04140479 0.2976573
# [3,]  0.02463362 0.7769173 0.8613714 0.94250880 0.8135982