library(ggplot2)

set.seed(1)
source("R/doPCA.R")
source("R/estimate_DR_curve.R")
source("R/kernel_smooth.R")
source("R/gcomp.R")
source("R/csPCA.R")

# Data generating processes -----------------------------------------------
sim_nonlinear_exp <- function(n=200, p=2, q=1, 
                           A=matrix(c(2,1),nrow=2), 
                           mu_C=1/2, Sigma_C=1, 
                           Sigma_e=diag(2), 
                           beta=c(1,.1), gamma=5, tau=1){
  
  p <- nrow(A); q <- ncol(A)
  e <- mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma_e)
  C <- mvrnorm(n = n, mu = mu_C, Sigma = Sigma_C)
  X <- 4*tanh(C %*% t(A) + e) 
  
  # Sigma_X <- A %*% Sigma_C %*% t(A) + Sigma_e
  
  omega <- beta/sqrt(sum(beta^2))
  Z     <- X %*% beta
  
  # Outcome model
  eps <- rnorm(n, mean = 0, sd = tau)
  Y   <- exp(-Z) + log(abs(C)) %*% gamma + eps
  
  return(list(Y = Y, X = X, Z = Z, C = C, 
              Sigma_C = Sigma_C, Sigma_e = Sigma_e, #Sigma_X = Sigma_X, 
              beta = beta, omega = omega, gamma = gamma, tau = tau))
}

sim_nonlinear <- function(n=100, p=2, q=1, 
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
  
  Y <- Z + C %*% gamma + eps
  
  return(list(Y = Y, X = X, Z=Z, C = C, 
              Sigma_C = Sigma_C, Sigma_e = Sigma_e, Sigma_X = Sigma_X, 
              beta = beta, omega = omega, gamma = gamma, tau = tau))
}

sim_linear <- function(n=100, p=2, q=1, 
                       A=matrix(c(2,1),nrow=2), 
                       mu_C=1/2, Sigma_C=1, 
                       Sigma_e=diag(c(9,1)), 
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
  
  return(list(Y = Y, X = X, C = C, 
              Sigma_C = Sigma_C, Sigma_e = Sigma_e, Sigma_X = Sigma_X, 
              beta = beta, omega = omega, gamma = gamma, tau = tau))
}


sim <- sim_nonlinear()

# Simulation --------------------------------------------------------------
B <- 100
list <- list()
for(i in 1:B){
  sim <- sim_linear()
  
  cs <- csPCA(sim$Y, sim$X, sim$C, maxit=1000)
  list[[i]] <- cs$omega
}

df <- do.call(rbind, list)
df[df[,1] < 0] <- -df[df[,1] < 0] 
colnames(df) <- c("x", "y")

domega <- sim$omega

# Add a label column to indicate which is the true vector
df$label <- "candidate"
df[nrow(df) + 1, ] <- c(domega[1], domega[2], "true")

# Convert label to factor so ggplot handles it correctly
df$label <- factor(df$label, levels = c("candidate", "true"))

# Plot vectors from origin (0, 0)
ggplot(df, aes(x = 0, y = 0, xend = x, yend = y)) +
  geom_segment(arrow = arrow(length = unit(0.2, "cm")), alpha = 0.3) +
  coord_equal() +
  theme_minimal() +
  labs(title = "Estimated Vectors and True Vector",
       x = "X", y = "Y") +
  theme(legend.title = element_blank())



# Likelihood evaluation ---------------------------------------------------

# Likelihood
f.do <- function(omega, Y, X, C, scaled = FALSE){
  print(omega)
  Z <- as.numeric( X %*% omega )
  DR <- gcomp(Y, Z, C)
  if(scaled){
    var(DR) / var(Z)
  } else{
    var(DR)
  }
}

mu_X <- gcomp(sim$Y, sim$X, sim$C)
f.cs <- function(omega, Y, X, C){
  print(omega)
  Z      <- as.numeric(X %*% omega)
    
  # E(mu(X) | Z = z) by kernel smoothing from Kennedy et al. (2017)
  reg <- kernel_smooth(mu_X, Z)$est
    
  # Optimizers are minimizers so take the negative of the variance
  var(reg)
}

# Evaluate at each theta
theta      <- seq(-pi, pi, length.out = 100)
omega_grid <- cbind(cos(theta), sin(theta))

# Calculate Objective for each PCA method
vals.do.scaled <- apply(omega_grid, 1, f.do, sim$Y, sim$X, sim$C, scaled = TRUE)
vals.do        <- apply(omega_grid, 1, f.do, sim$Y, sim$X, sim$C)
vals.cs        <- apply(omega_grid, 1, f.cs, sim$Y, sim$X, sim$C)

# omega.hat.do <- omega_grid[which.max(vals.do),] |> get_theta()
omega.do        <- optim_dopca(sim$beta, sim$Sigma_e, sim$Sigma_X, scaled = FALSE)
omega.do.scaled <- optim_dopca(sim$beta, sim$Sigma_e, sim$Sigma_X, scaled = TRUE)


dt <- data.frame(theta = theta, doPCA = vals.do, csPCA = vals.cs, doPCA_scaled = vals.do.scaled) |>
  tidyr::pivot_longer(
    cols = c(doPCA, doPCA_scaled, csPCA),
    names_to = "method"
  )

th_true <- get_theta(sim$omega)   # numeric x-location for the vline
y_top   <- max(dt$value, na.rm = TRUE)  # or choose a nicer offset

ggplot(dt, aes(theta, value)) +
  geom_line(aes(group = method, colour = method)) +
  geom_vline(xintercept = th_true, linetype = "dashed") +
  geom_vline(xintercept = get_theta(omega.do), linetype = "dashed") +
  geom_vline(xintercept = get_theta(omega.do.scaled), linetype = "dashed") +
  annotate("text",
           x = th_true,
           y = y_top,
           label = "true~omega",  # plotmath expression string
           parse = TRUE,
           vjust = -0.3,          # nudge above the top point; adjust as needed
           colour = "red") +
  ggtitle("Sample objective function for n=100, linear, Sigma_e = diag(9,1)") +
  xlab(expression(theta)) +
  theme_bw()

ggsave("results/objective_vis_linear_high_varX1.png")


obj <- list(theta = theta, 
            omega_grid = omega_grid, 
            vals.do.scaled = vals.do.scaled,
            vals.do = vals.do,
            vals.cs = vals.cs
      )

saveRDS(obj, "sim_nonlinear_exp_obj.rds")
