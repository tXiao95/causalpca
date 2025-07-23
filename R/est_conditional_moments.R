est_conditional_moments <- function(Y, X, C) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  C <- as.data.frame(C)
  
  # Build design matrix with intercept
  Cmm <- model.matrix(~ . + 0, data = C)      # n x k
  qrC <- qr(Cmm)
  n   <- nrow(Cmm); k <- qrC$rank
  
  # Residuals (handles multivariate RHS)
  EX <- qr.resid(qrC, X)                  # n x p
  EY <- qr.resid(qrC, Y)                  # n x q
  
  # Cross-products
  Var_X_given_C  <- crossprod(EX) / (n - k)            # p x p
  Cov_XY_given_C <- crossprod(EX, EY) / (n - k)        # p x q
  
  list(
    Var_X_given_C  = Var_X_given_C,
    Cov_XY_given_C = Cov_XY_given_C,
    residuals_X = EX,
    residuals_Y = EY,
    df = n - k,
    design = Cmm
  )
}
