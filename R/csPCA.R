source("R/misc.R")
source("R/kernel_smooth.R")
source("R/gcomp.R")

library(MAVE)

csPCA <- function(Y, X, C, maxit = 5000, verbose = FALSE, omega0=NULL){
  n <- nrow(X); p <- ncol(X); q <- ncol(C)
  
  # ---------------------------------------------------------------------------
  # Objective as a function of unconstrained theta
  # ---------------------------------------------------------------------------
  
  # Estimate mu(X) for multidimensional continuous treatment (g-comp)
  mu_X <- gcomp(Y, X, C)
  
  objective_omega <- function(omega){
    # We enforce the unit norm in the optimization here
    omega  <- omega / sqrt(sum(omega^2))
    Z      <- as.numeric(X %*% omega)
    
    # E(mu(X) | Z = z) by kernel smoothing, borrowed from Kennedy et al. (2017)
    reg <- kernel_smooth(mu_X, Z)$est
    
    # Optimizers are minimizers so take the negative of the variance
    val    <- -var(reg)
    
    cat('Omega:', omega, '\n')
    cat('Norm: ', sqrt(sum(omega^2)), '\n')
    cat('Objective:', val, '\n\n')
    
    return(val)
  }
  
  constraint_omega_eq <- function(omega){
    constr <- sum(omega^2) - 1
    jac    <- 2 * omega
    list("constraints" = constr,
         "jacobian"    = jac)
  }
  
  # ---------------------------------------------------------------------------
  # Optimization routine
  # ---------------------------------------------------------------------------
  # initial value: first PC of X
  if(is.null(omega0)){
    omega0 <- prcomp(X, center = TRUE, scale. = FALSE)$rotation[,1]
  }
  opts   <- list(algorithm = "NLOPT_LN_NELDERMEAD",
                maxeval = maxit,
                xtol_rel = 1e-4,
                print_level = 1)
  
  # Solve optimization
  opt <- nloptr(
    x0 = omega0,
    eval_f = objective_omega,
    # eval_g_eq = constraint_omega_eq,
    # Bounds on each element of the omega. Since unit norm, these must be true. 
    lb        = rep(-1, length(omega0)),
    ub        = rep( 1, length(omega0)),
    opts = opts
  )
  
  # Get solution
  omega_opt <- opt$solution |> unitvec()
  Z_opt     <- as.numeric(X %*% omega_opt)
  
  # Return solution, dose-response curve, and optimization metadata
  res <- list(omega       = omega_opt,
              Z           = Z_opt,
              opt         = opt)
  
  class(res) <- 'csPCA'
  return(res)
}
