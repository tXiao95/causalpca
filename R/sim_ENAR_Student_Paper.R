library(MASS)

# --- helpers ---
softplus <- function(x) log1p(exp(x))  # smooth >=0
# (optional) target a given positive mean for log-normal: mu_L = log(m) - 0.5*s2

simulate_data_pos <- function(n, p = 20, q,
                              rho = 0.1,
                              lambda = 0.0,
                              Theta,                 # q x (p-4)
                              beta,                  # 20 x 2 from make_beta(...)
                              alpha = rep(1, q),
                              method = c("softplus", "lognormal", "gamma")) {
  method <- match.arg(method)
  
  # Covariates
  C <- MASS::mvrnorm(n = n, mu = rep(0, q), Sigma = diag(q))
  m1_C <- rowSums(C)
  m2_C <- as.vector(C %*% ((-1)^(1:q)))
  
  # --- Latent Gaussian versions of X1..X4 (keep your structure) ---
  X1_lat <- rnorm(n, m1_C, sd = 1)
  X2_lat <- rnorm(n, m2_C, sd = 1)
  X3_lat <- rnorm(n, mean = abs(X1_lat + X2_lat),
                  sd   = 0.1)
  X4_lat <- rnorm(n, mean = sqrt(abs(X1_lat + X2_lat) + 0.1),
                  sd   = 0.1)
  
  # --- Last 16 with AR(1) dependence on a latent Gaussian, mean linked to C ---
  Sigma_X <- toeplitz(rho^(0:((p-4)-1)))
  mu_lat  <- t(Theta) %*% t(C)            # (p-4) x n
  X5_20_lat <- t(apply(mu_lat, 2, function(mu_i)
    MASS::mvrnorm(1, mu = mu_i, Sigma = Sigma_X)
  ))  # n x (p-4)
  
  # --- Map to positive reals ---
  make_pos <- switch(method,
                     "softplus" = function(z) softplus(z),   # smoother, lighter tail
                     "lognormal" = function(z) exp(z),       # heavier tail, strictly > 0
                     "gamma" = function(mu) {                # Gamma with mean mu and fixed cv
                       # Choose coefficient of variation (cv); smaller cv -> less noise
                       cv <- 0.5
                       shape <- 1/(cv^2)
                       scale <- (mu / shape)
                       rgamma(length(mu), shape = shape, scale = scale)
                     }
  )
  
  # For gamma, we want a positive mean to feed; take positive transform of latents first
  if (method == "gamma") {
    X1 <- make_pos(softplus(X1_lat))
    X2 <- make_pos(softplus(X2_lat))
    X3 <- make_pos(softplus(X3_lat))
    X4 <- make_pos(softplus(X4_lat))
    X5_20 <- apply(X5_20_lat, 2, function(col) make_pos(softplus(col)))
    X5_20 <- matrix(X5_20, nrow = n, ncol = p - 4)
  } else {
    X1 <- make_pos(X1_lat)
    X2 <- make_pos(X2_lat)
    X3 <- make_pos(X3_lat)
    X4 <- make_pos(X4_lat)
    X5_20 <- make_pos(X5_20_lat)
  }
  
  X <- cbind(X1, X2, X3, X4, X5_20)
  
  # --- Outcome ---
  Z <- X %*% beta                   # n x 2
  g_C <- as.vector( abs(C %*% alpha) + 0.3 * C[, 1]^2 - 0.3 * sin(C[, 2]))
  h_Z <- Z[, 1] + Z[, 2]
  
  # Default: unconstrained Y
  Y <- (h_Z + g_C + lambda * rowSums(X) * rowSums(C) + rnorm(n, 0, 1) )
  
  list(C = C, X = X, Y = Y)
}

data <- simulate_data_pos(n = 100, p = 20, q = 10, 
                          rho = 0.6, lambda = 0.01, 
                          Theta = make_Theta(q = 10), beta = make_beta())

# Show correlation between exposures for simulation. 
# Compare to correlation in the data. 
corrplot::corrplot(cor(data$X))

# True dose-response is E(Y(x)) = h(Z) = Z[,1] + Z[, 2] = beta1'X + beta2'X
# Report true dose-response curve once decided (just pick one case for stuent paper). 
