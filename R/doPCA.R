library(nloptr)

source("R/estimate_DR_curve.R")
source("R/misc.R")

#' Causal doPCA
#'
#' @description Estimate the first causal principal direction using doPCA
#'
#' @param Y Numeric vector (n x 1) – outcome.
#' @param X Numeric matrix (n x p) – multidimensional treatment.
#' @param C Numeric matrix (n x q) – baseline confounders.
#' @param maxit Integer.  Maximum iterations 
#' @param verbose Logical.  Print objective value during optimisation.
#'
#' @return List with elements
#'
#' @export
#'
#' @examples
#' set.seed(1)
# n <- 200; p <- 3; q <- 2
# X <- matrix(rnorm(n*p), n, p)
# C <- matrix(rnorm(n*q), n, q)
# omega_true <- c(1,0.5,-0.5); omega_true <- omega_true/sqrt(sum(omega_true^2))
# Z <- X %*% omega_true
# Y <- 2*Z + 0.5*C[,1] + rnorm(n)
# res <- doPCA(Y, X, C)
# res$omega
doPCA <- function(Y, X, C, mu_est = "gcomp", scaled = FALSE, maxit = 5000, verbose = FALSE, omega0 = NULL){
  n <- nrow(X); p <- ncol(X); q <- ncol(C)
  
  # ---------------------------------------------------------------------------
  # Objective as a function of unconstrained theta
  # ---------------------------------------------------------------------------
  objective_omega <- function(omega){
    # We enforce the unit norm in the optimization here.
    omega  <- omega / sqrt(sum(omega^2))
    Z      <- as.numeric(X %*% omega)
    if(mu_est == "gcomp"){
      mu_hat <- gcomp(Y, Z, C)
    }
    if(mu_est == "DR"){
      mu_hat <- estimate_DR_curve(Y, Z, C, Z.new = Z)$DR_curve
    }
    # Optimizers are minimizers so take the negative of the variance
    val    <- -var(mu_hat)
    if(scaled){val <- val / var(Z)}
    
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
    # return(constr)
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
    # Bounds on each element of the omega
    lb        = rep(-1, length(omega0)),
    ub        = rep( 1, length(omega0)),
    opts = opts
  )
  
  # Get solution
  omega_opt <- opt$solution |> unitvec()
  Z_opt     <- as.numeric(X %*% omega_opt)
  # mu_opt    <- estimate_DR_curve(Z_opt, Y, C, Z.new = Z_opt)
  
  # Return solution, dose-response curve, and optimization metadata
  res <- list(omega       = omega_opt,
              # mu_hat      = mu_opt,
              Z           = Z_opt,
              opt = opt)
  
  class(res) <- 'doPCA'
  return(res)
}