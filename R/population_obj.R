library(nloptr)

optim_dopca <- function(beta, Sigma_e, Sigma_X, scaled = FALSE) {
  ## ---- analytic gradient of the objective -------------------
  f <- function(w) {
    a <- as.numeric(beta %*% Sigma_e %*% w)       # β'Σe ω   (scalar)
    b <- as.numeric(w %*% Sigma_e %*% w)          # ω'Σe ω   (scalar)
    
    if (scaled) {
      val  <- (a^2)/(b^2)
    } else {
      wXw <- w %*% Sigma_X %*% w         # ω'Σxx ω
      f1  <- (a^2)/(b^2)                                      # inner ratio
      val <- f1 * wXw
    }
    return(-val)
  }
  
  ## ---- equality constraint  g(ω) = ‖ω‖² - 1 = 0 ------------
  geq <- function(w) sum(w^2) - 1
  
  ## ---- initial guess: unit vector along β ------------------
  w0 <- beta / sqrt(sum(beta^2))
  
  sol <- nloptr(
    x0        = w0,
    eval_f    = f,        
    eval_g_eq = geq,
    opts      = list(
      algorithm = "NLOPT_LN_COBYLA",
      xtol_rel  = 1e-8,
      maxeval   = 1000
    )
  )
  
  sol$solution / sqrt(sum(sol$solution^2))     # ensure exact unit norm
}

optim_cspca <- function(beta, Sigma_XX) {
  ## ---- analytic gradient of the objective -------------------
  feval <- function(w) {
    a <- as.numeric(beta %*% Sigma_XX %*% w)          # β'Σe ω   (scalar)
    b <- as.numeric(w %*% Sigma_XX %*% w)          # ω'Σe ω   (scalar)
      
    f  <- t(beta) %*% Sigma_XX %*% beta - a^2 / b
    return(f)
  }
  
  ## ---- equality constraint  g(ω) = ‖ω‖² - 1 = 0 ------------
  geq <- function(w) sum(w^2) - 1
  
  ## ---- initial guess: unit vector along β ------------------
  w0 <- beta / sqrt(sum(beta^2))
  
  sol <- nloptr(
    x0        = w0,
    eval_f    = feval,        # returns both objective and gradient
    eval_g_eq = geq,
    opts      = list(
      algorithm = "NLOPT_LN_COBYLA",
      xtol_rel  = 1e-8,
      maxeval   = 1000
    )
  )
  
  # return(sol)
  sol$solution / sqrt(sum(sol$solution^2))     # ensure exact unit norm
}