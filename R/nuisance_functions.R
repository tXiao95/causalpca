library(SuperLearner)

# Outcome regression and gcomp --------------------------------------------

fit_outcome_regression <- function(Y, V, ...) {
  stopifnot(length(Y) == nrow(V))
  
  V <- as.data.frame(V)
  p <- ncol(V)
  
  sl_fit <- SuperLearner::SuperLearner(Y = Y, X = V, ...)
  
  out <- list(
    sl_fit = sl_fit,
    p = p,
    V = V            # store predictors actually used
  )
  
  class(out) <- "or_fit"
  out
}

gcomp_from_fit <- function(or_fit, X.new = NULL) {
  stopifnot(inherits(or_fit, "or_fit"))
  sl_fit <- or_fit$sl_fit
  C <- or_fit$C
  n <- nrow(C)
  p <- or_fit$p
  
  # default: evaluate at observed X
  if (is.null(X.new)) {
    X.new <- or_fit$X
  } else {
    X.new <- as.data.frame(X.new)
    if (ncol(X.new) != p) stop("X.new must have ", p, " columns.")
  }
  
  m <- nrow(X.new)
  
  gcomp_est <- vapply(seq_len(m), function(i) {
    Xi <- X.new[i, , drop = FALSE]
    df_new <- cbind(Xi[rep(1, n), , drop = FALSE], C)
    
    # Depending on SL version, one of these will work:
    # pred <- SuperLearner::predict.SuperLearner(sl_fit, newdata = df_new)$pred
    pred <- predict(sl_fit, newdata = df_new)$pred
    
    mean(pred)
  }, numeric(1L))
  
  names(gcomp_est) <- paste0("x", seq_len(m))
  gcomp_est
}

# convenience wrapper (fit if not provided)
gcomp <- function(Y, X, C, X.new = NULL, ..., or_fit = NULL) {
  if (is.null(or_fit)) {
    or_fit <- fit_outcome_regression(Y, X, C, ...)
  }
  gcomp_from_fit(or_fit, X.new = X.new)
}

# Nuisances (no beta) -----------------------------------------------------
SL.lib <- c("SL.glm", "SL.ranger")

or_fit <- fit_outcome_regression(Y, X, C, SL.library = SL.lib)
muX    <- gcomp_from_fit(or_fit)   # X.new defaults to original X
varX   <- gcomp(Y^2, X, C, SL.library = SL.lib) - muX


# Data --------------------------------------------------------------------
compute_score_eq7 <- function(
    Y,                   # n 
    X,                   # n x p
    C,                   # n x q
    EY_XC,               # n (outcome regression)
    pX,                  # n (stabilisation weights)
    pXC,                 # n (propensity score)
    m_hat,               # n : m̂(β^T X_i)
    mprime_hat,          # n x d : m̂'(β^T X_i)
    sigma2_hat,          # n : σ̂^2(X_i)
    EX_over_sigma2,      # n x p : Ê{X/σ̂^2(X) | β^T X_i}
    E_inv_sigma2         # n : Ê{1/σ̂^2(X) | β^T X_i}
) {
  stopifnot(is.matrix(X))
  n <- nrow(X); p <- ncol(X)
  stopifnot(length(Y) == n,
            length(m_hat) == n,
            nrow(mprime_hat) == n,
            length(sigma2_hat) == n,
            nrow(EX_over_sigma2) == n,
            ncol(EX_over_sigma2) == p,
            length(E_inv_sigma2) == n)
  
  d <- ncol(mprime_hat)
  
  # Propensity weight that is n x 1 vector
  W <- pX / pXC
  
  # p-vector inside brackets for each i
  # g_i = (1/sigma2_hat_i) * bracket_i
  ratio_mat <- EX_over_sigma2 / E_inv_sigma2  # n x p (recycles E_inv_sigma2 by row)
  bracket   <- X - ratio_mat                  # n x p 
  g_mat     <- bracket / sigma2_hat               # n x p (recycles sigma2_hat by row, each row gets same sigma)
  
  # residuals
  r_Y  <- Y - m_hat                              
  r_EY <- EY_XC - m_hat
  
  # Pass 1: compute vec(Yterm) and vec(EYterm) which are n x (pd) matrices
  Yterm_vec  <- matrix(0, nrow = n, ncol = p*d)
  EYterm_vec <- matrix(0, nrow = n, ncol = p*d)
  
  for (i in 1:n) {
    temp            <- (g_mat[i, ] %o% mprime_hat[i, ]) * W[i]  # p x d
    Yterm_vec[i, ]  <- as.vector(temp * r_Y[i])      # vec(p x d)
    EYterm_vec[i, ] <- as.vector(temp * r_EY[i])     # vec(p x d)
  }
  
  # Fit E(phi|C) ≈ E(EYterm | C) componentwise, also n x (pd)
  Ephi_vec <- fit_pd_list(Ymat = EYterm_vec, V = C)
  
  # Compute score for each i, its sample average, and the cross product. 
  S_i_vec    <- Yterm_vec - EYterm_vec + Ephi_vec
  S_bar_vec  <- colMeans(S_i_vec)
  S2_bar     <- crossprod(S_i_vec) / n
  
  return(list(S_bar_vec = S_bar_vec, S2_bar = S2_bar))
}

fit_pd_list <- function(Ymat, V, ...) {
  Y <- as.matrix(Ymat)
  V <- as.matrix(V)
  stopifnot(nrow(Y) == nrow(V))
  
  out <- vector("list", ncol(Y))
  
  for (j in seq_len(ncol(Y))) {
    fit <- fit_outcome_regression(Y = Y[,j], V = V, ...)
    out[[j]] <- ( predict( fit$sl_fit) )$pred
  }
  
  do.call(cbind, out)
}

# Test --------------------------------------------------------------------

set.seed(123)

# Dimensions
n <- 100
p <- 10
q <- 4
d <- 2

# Basic data
Y <- rnorm(n)

X <- matrix(rnorm(n * p), n, p)
C <- matrix(rnorm(n * q), n, q)

# Nuisance-type quantities
EY_XC <- rnorm(n)              # outcome regression prediction

pX  <- runif(n, 0.2, 0.8)      # stabilisation weights
pXC <- runif(n, 0.2, 0.8)      # propensity scores

m_hat      <- rnorm(n)                 # m̂(β^T X_i)
mprime_hat <- matrix(rnorm(n * d), n, d)

sigma2_hat <- abs(rnorm(n)) + 0.5      # ensure positive

EX_over_sigma2 <- matrix(rnorm(n * p), n, p)
E_inv_sigma2   <- abs(rnorm(n)) + 0.5  # positive

# Now test the function
out <- compute_score_eq7(
  Y = Y,
  X = X,
  C = C,
  EY_XC = EY_XC,
  pX = pX,
  pXC = pXC,
  m_hat = m_hat,
  mprime_hat = mprime_hat,
  sigma2_hat = sigma2_hat,
  EX_over_sigma2 = EX_over_sigma2,
  E_inv_sigma2 = E_inv_sigma2
)

str(out)

MASS::ginv(out$S2_bar) %*% out$S_bar_vec

