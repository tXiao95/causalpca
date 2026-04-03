simulate_data_ENAR <- function(n,
                          p = 20, 
                          q = 10,
                          rho = 0.8,
                          Z1_coef = 1,
                          Z2_coef = 0.5,
                          Z12_coef = 0.5,
                          h_Z_coef = 0.1,
                          g_C_coef = 0.1,
                          var_scale = 3,
                          interaction_coef = 0.1) {
  if(p < 6){
    stop("Number of exposure variables must be greater than or equal to 6")
  }
  # Confounders [C] ---------------------------------------------------------
  Sigma_C <- toeplitz(c(1, .5, rep(0, q-2)))
  C       <- MASS::mvrnorm(n = n, mu = rep(0, q), Sigma = Sigma_C)

  # Inflate variance in X's that don't feed into beta to distract PCA 
  idx_beta  <- c(1:2, 7:10)
  scale_vec <- rep(1, p)
  scale_vec[-idx_beta] <- var_scale
  D         <- diag(scale_vec, p, p)
  
  # Exposures [X | C] -------------------------------------------------------
  Sigma_X <- D %*% toeplitz(rho^(0:((p)-1))) %*% D
  Theta   <- matrix(  ( 1:p / (1:q)^2 ), nrow = q, ncol = p)
  X_mean  <- pnorm(C %*% Theta) 
  #X_mean  <- (C %*% Theta) 
  
  X <- t(apply(X_mean, 1, function(mu_i)
    MASS::mvrnorm(1, mu = mu_i, Sigma = Sigma_X)
  ))  # n x (p-4)
  
  # 2 dimensions, first direction is 4, second is 2. 
  beta <- matrix(c(rep(1 / sqrt(2), 2),rep(0, p - 2), 
                 c(rep(0, p - 4), rep(1 / sqrt(4), 4)) ),nrow = p, ncol = 2)
  # the unique Projection matrix from R^p to R^d
  P_beta <- beta %*% solve(t(beta) %*% beta) %*% t(beta)
  
  # Dimension reduction
  Z   <- X %*% beta                   # n x 2
  
  # --- Outcome ---
  g_C   <- sin( C %*% Theta[,1] ) 
  h_Z   <- Z1_coef * Z[,1] + Z2_coef * Z[,2]^2 + Z12_coef * Z[,1]*Z[,2]
  eps_Y <- rnorm(n, mean = 0, sd = 0.5)
  
  # Y response model
  Y     <- h_Z_coef * h_Z + g_C_coef * g_C + interaction_coef * (Z[,1]*C[,1]) + eps_Y * (sqrt( 0.5 + pnorm(rowSums(C[,1:2])) )) 
  
  # Causal mean at observed points (g_C is mean zero)
  mu_X <- h_Z_coef * h_Z
  
  # Causal mean function
  h_Z_fun <- function(Z1, Z2){
    return( h_Z_coef * (Z1_coef * Z1 + Z2_coef * Z2^2 + Z12_coef * Z1 * Z2 ) )
  }
  
  list(Y = Y, C = C, X = X, mu_X = mu_X, beta_true = beta, P_beta= P_beta, h_Z_fun = h_Z_fun, d = ncol(beta))
}

#' Minimal First-Principles Causal SDR Simulation (With Correlation & Interactions)
#' 
#' @param n Sample size
#' @param p Number of exposures (default 10)
#' @param q Number of confounders (default 5)
#' @param noise_sd Noise level for the outcome
#' @param rho_X Autoregressive correlation between X variables
#' @param interaction_coef Strength of the interaction (Set to 0 for RP to work, >0 for RP to fail)
#' @param weak_dim_signal whether to use the simulation that results in a weak dimension signal
#' @param sparse whether to use the sparse assumption (only select X1 and X2)
#' @return A list of data and true structural components

# Simulate strong dimension signal
simulate_causal_sdr_simple <- function(n = 1000, 
                                       p = 10, 
                                       q = 5, 
                                       noise_sd = 0.5,
                                       rho_X = 0.8,
                                       interaction_coef = 5.0, 
                                       weak_dim_signal = FALSE,
                                       sparse = TRUE) {
  
  if (p < 6) stop("p must be >= 6 to accommodate causal, MAVE traps, and PCA traps.")
  if (q < 5) stop("q must be >= 5 to accommodate the confounding structures and interactions.")
  
  # 1. Confounders (Independent Standard Normal)
  C <- matrix(rnorm(n * q), nrow = n, ncol = q)
  colnames(C) <- paste0("C", 1:q)
  
  # 2. Exposures (X) with AR(1) Correlation
  R <- toeplitz(rho_X^(0:(p-1))) # Base AR(1) Correlation Matrix
  
  # Scale variance: Massive variance for X5 and X6 to bait PCA
  scale_vec <- rep(1, p)
  scale_vec[5:6] <- 10.0 
  D <- diag(scale_vec, p, p)
  
  Sigma_X <- D %*% R %*% D # Final Covariance Matrix
  
  # Mean of X given C
  mu_X_cond <- matrix(0, nrow = n, ncol = p)
  mu_X_cond[, 1] <- 0.5 * C[, 1]
  mu_X_cond[, 2] <- 0.5 * C[, 2]
  mu_X_cond[, 3] <- 3.0 * C[, 3] # MAVE Trap 1
  mu_X_cond[, 4] <- 3.0 * C[, 4] # MAVE Trap 2
  # X5 through Xp have conditional mean 0
  
  X <- mu_X_cond + MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma_X)
  colnames(X) <- paste0("X", 1:p)
  
  # 3. True Causal Directions (beta)
  # The true causal subspace spans exactly X1 and X2. (d = 2)
  beta <- matrix(0, nrow = p, ncol = 2)
  if(weak_dim_signal){
    beta[1:2, 1] <- 1 / sqrt(2)
    beta[3, 2] <- 1 
  } else if(sparse){
    beta[1,1] <- 1
    beta[2,2] <- 1
  } else{
    # Z1 uses odd-indexed variables: X1, X3, X7, X9
    beta[c(1, 3, 7, 9), 1] <- 0.5  # which is exactly 1/sqrt(4)
    
    # Z2 uses even-indexed variables: X2, X4, X8, X10
    beta[c(2, 4, 8, 10), 2] <- 0.5 # which is exactly 1/sqrt(4)
  }
  
  Z <- X %*% beta  
  
  # 4. Causal Exposure-Response Surface: h(Z)
  h_Z_fun <- function(Z1, Z2) {
    #return( 4 * tanh(Z1) + 2 * Z2^2 + Z1*Z2) # THIS IS FOR THE NNET_EASY STUDY
    return( 4 * sin(Z1) + 2 * Z2^2 + Z1*Z2)  #THIS IS FOR THE MAIN STUDY WITH X1 and X2 (ENAR)
    #return( 3 * tanh(Z1) + 8 * (pnorm(2 * Z2) - 0.5) + 0.2 * Z2*Z1 )
  }
  mu_causal <- h_Z_fun(Z[, 1], Z[, 2])
  
  # # 5. Additive Confounding Effect: g(C)
  g_C <- 5.0 * tanh(C[, 1]) + 
         5.0 * (C[, 2]^2 - 1) + 
         5.0 * (C[, 3]^2 - 1) + 
         #2.5 * C[, 4]             # THIS IS FOR THE NNET_EASY STUDY
         5.0 * sin(1.5 * C[, 4]) #THIS IS FOR THE MAIN STUDY WITH X1 AND X2 (ENAR)
  
  # 6. Interaction Trap (Breaks RP when > 0)
  interaction_term <- interaction_coef * C[, 5] * (Z[, 1] + Z[, 2])
  
  # 7. Final Outcome Y
  Y <- mu_causal + g_C + interaction_term + rnorm(n, 0, sd = noise_sd)
  
  # True Projection Matrix
  P_beta <- beta %*% solve(t(beta) %*% beta) %*% t(beta)
  
  return(list(
    Y = Y, C = C, X = X, Z = Z, 
    mu_X = mu_causal, 
    beta_true = beta, 
    P_beta = P_beta, 
    h_Z_fun = h_Z_fun, 
    d = ncol(beta)
  ))
}