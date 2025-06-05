# A simple illustration to show the deficiences of PCA, SPCA, CCA, and their partial counterparts
# for causal dimension reduction.

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

# ========================
# Example: causal setting with p = 3 and confounding
# ========================

set.seed(42)
n <- 200; p <- 3; q <- 1  # 10 exposures, 1 confounder

# Step 1: Generate confounder(s)
C <- matrix(rnorm(n * q), n, q)

# Step 2: Generate exposures X with some depending on C
X     <- matrix(rnorm(n * p), n, p)
X[,1] <- 0.2 * C[,1] + rnorm(n)     # Mild confounding + signal
X[,2] <- 5 * C[,1] + rnorm(n)    # Strong confounding + no signal
                                # X3 is just independent noise

# Step 3: Generate outcome Y as a function of X1 and C. 
Y_linear   <- ( X[,1] + C[,1] + rnorm(n) )
Y_threshold <- X[,1] * ifelse(C[,1] > 0, 1, -1) + rnorm(n)

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
# PCA        SPCA         CCA      P-SPCA         P-CCA
# [1,] -0.04325582  0.74976672 0.767785922  0.99849388  0.9691108604
# [2,] -0.99898813  0.66154747 0.640706359 -0.02987290 -0.2466249555
# [3,]  0.01231462 -0.01431151 0.000373572  0.04601715 -0.0005211904
# 
# $`Threshold/nonlinear confounding\n                            `
# PCA      SPCA       CCA     P-SPCA     P-CCA
# [1,] -0.04325582 0.2776127 0.1547308 0.21095242 0.1125317
# [2,] -0.99898813 0.4938294 0.5146281 0.05236491 0.8367025
# [3,]  0.01231462 0.8240532 0.8433364 0.97609272 0.5359716
# 
