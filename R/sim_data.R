library(MASS)

sim_unif_e <- function(n=100, p=2, q=1, 
                           A=matrix(c(2,1),nrow=2), 
                           mu_C=1/2, Sigma_C=1, 
                           Sigma_e=diag(2), 
                           beta=c(1,.1), gamma=5, tau=1){
  
  p <- nrow(A); q <- ncol(A)
  e <- matrix(runif(n*p, -3, 3), nrow = n)
  C <- mvrnorm(n = n, mu = mu_C, Sigma = Sigma_C)
  X <- C %*% t(A) + e
  
  Sigma_X <- A %*% Sigma_C %*% t(A) + Sigma_e
  
  omega <- beta/sqrt(sum(beta^2))
  Z          <- X %*% beta
  
  eps <- rnorm(n, mean = 0, sd = tau)
  
  Y <- Z^2 + C %*% gamma + eps
  
  return(list(Y = Y, X = X, Z=Z, C = C, 
              Sigma_C = Sigma_C, Sigma_e = Sigma_e, Sigma_X = Sigma_X, 
              beta = beta, omega = omega, gamma = gamma, tau = tau))
}

sim_linear <- function(n=100, p=2, q=1, 
                       A=matrix(c(2,1),nrow=2), 
                       mu_C=1/2, Sigma_C=1, 
                       Sigma_e=diag(c(1,1)), 
                       beta=c(1,.1), gamma=5, tau=1){
  
  p <- nrow(A); q <- ncol(A)
  e <- mvrnorm(n = n, mu = rep(0, p), Sigma_e)
  C <- mvrnorm(n = n, mu = mu_C, Sigma = Sigma_C)
  X <- C %*% t(A) + e
  
  Sigma_X <- A %*% Sigma_C %*% t(A) + Sigma_e
  
  omega <- beta/sqrt(sum(beta^2))
  Z          <- X %*% beta
  
  eps <- rnorm(n, mean = 0, sd = tau)
  
  Y <- Z + C %*% gamma + eps
  
  return(list(Y = Y, X = X, Z = Z, C = C, 
              Sigma_C = Sigma_C, Sigma_e = Sigma_e, Sigma_X = Sigma_X, 
              beta = beta, omega = omega, gamma = gamma, tau = tau))
}