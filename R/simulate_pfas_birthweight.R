library(MASS)

simulate_pfas_birthweight <- function(n = 303, 
                                      #noise_sd = 245, 
                                      noise_sd = 200, 
                                      #interaction_coef = 10.0) {
                                      interaction_coef = 200) {
  
  p <- 4 # PFOS, PFOA, PFNA, PFHxS
  q <- 6 # Confounders
  
  # ---------------------------------------------------------
  # 1. Confounders (C)
  # ---------------------------------------------------------
  C <- matrix(rnorm(n * q), nrow = n, ncol = q)
  colnames(C) <- paste0("C", 1:q)
  
  # ---------------------------------------------------------
  # 2. PFAS Exposures (X) with Confounding
  # ---------------------------------------------------------
  # Empirical correlation matrix from the real analysis dataset
  Cor_X <- matrix(c(
    1.0000000, 0.5800987, 0.4456830, 0.3531313,
    0.5800987, 1.0000000, 0.6703732, 0.3366613,
    0.4456830, 0.6703732, 1.0000000, 0.4118724,
    0.3531313, 0.3366613, 0.4118724, 1.0000000
  ), nrow = 4, ncol = 4, byrow = TRUE)
  
  # Scale base variance down slightly (0.8) so that after adding 
  # the confounder variance below, the final Z stays mostly in [-2, 2]
  Sigma_X <- Cor_X * 0.8 
  
  # Conditional mean of X given C (This embeds the confounding!)
  mu_X_cond <- matrix(0, nrow = n, ncol = p)
  mu_X_cond[, 1] <- 0.4 * C[, 1] + 0.2 * C[, 2] # PFOS
  mu_X_cond[, 2] <- 0.4 * C[, 2] - 0.2 * C[, 3] # PFOA
  mu_X_cond[, 3] <- 0.3 * C[, 1] + 0.3 * C[, 4] # PFNA
  mu_X_cond[, 4] <- 0.5 * C[, 3]                # PFHxS
  
  # Generate X
  X <- mu_X_cond + mvrnorm(n, mu = rep(0, p), Sigma = Sigma_X)
  colnames(X) <- c("PFOS", "PFOA", "PFNA", "PFHxS")
  
  # ---------------------------------------------------------
  # 3. True Causal Directions (beta) and Latent Factor (Z)
  # ---------------------------------------------------------
  # Desired proportions for the true causal mechanism
  raw_beta <- c(0.6, 0.2, 0.2, -0.6)
  
  # Normalize to create a true orthonormal unit vector (L2 norm = 1)
  beta <- matrix(raw_beta / sqrt(sum(raw_beta^2)), ncol = 1) 
  
  # Calculate Z (d = 1)
  Z <- X %*% beta
  colnames(Z) <- "Z1"
  
  # ---------------------------------------------------------
  # 4. Causal Exposure-Response Surface: h(Z)
  # ---------------------------------------------------------
  h_Z_fun <- function(Z_val) {
    #return( 3300 - 80 * tanh(Z_val + 0.5) )
    return( 3300 - 150 * tanh(Z_val + 0.5) )
  }
  mu_causal <- h_Z_fun(Z[, 1])
  
  # ---------------------------------------------------------
  # 5. Additive Confounding Effect on Outcome: g(C)
  # ---------------------------------------------------------
  # g_C <- 80 * C[, 1] + 
  #        160 * (C[, 2]^2 - 1) + 
  #        160 * (C[, 3]^2 - 1) + 
  #        80 * C[, 4]
  g_C <- 50 * C[, 1] + 
    80 * pmax(C[, 2] - 0.5, 0) +  
    50 * abs(C[, 3]) +            
    50 * C[, 4]
  # ---------------------------------------------------------
  # 6. Interaction Trap
  # ---------------------------------------------------------
  interaction_term <- interaction_coef * C[, 5] * Z[, 1]
  
  # ---------------------------------------------------------
  # 7. Final Outcome (Y)
  # ---------------------------------------------------------
  Y <- mu_causal + g_C + interaction_term + rnorm(n, 0, sd = noise_sd)
  
  # True Projection Matrix
  P_beta <- beta %*% solve(t(beta) %*% beta) %*% t(beta)
  
  return(list(
    Y = Y, 
    C = C, 
    X = X, 
    Z = Z, 
    mu_X = mu_causal, 
    beta_true = beta, 
    P_beta = P_beta, 
    h_Z_fun = h_Z_fun, 
    d = 1
  ))
}
