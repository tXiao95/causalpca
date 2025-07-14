library(MASS)

# For B = 100 iterations
# 1. Solve doPCA
# 2. Record all tracked solutions and objective evaluations (from beginning to end)
# 3. Plot the likelihood
# 4. Compare to the true beta, doPCA truth, and sample estimation

# Do the same for p = 3

source("R/doPCA.R")
source("R/estimate_DR_curve.R")
source("R/population_obj.R")

set.seed(1)


# Settings
n <- 200; p <- 2; q <- 2

A <- matrix(rnorm(p * q, 1, 1), nrow = p, ncol = q)

# Confounders
Sigma_C <- diag(q)
mu_C    <- rep(1, q)
C       <- mvrnorm(n = n, mu = mu_C, Sigma_C)

# Treatment
Sigma_e <- diag(p)
mu_X    <- A %*% mu_C 
Sigma_X <- A %*% Sigma_C %*% t(A) + Sigma_e

X <- mvrnorm(n = n, mu_X, Sigma_X)

# True beta and Y
beta        <- c(1,0.5) 
beta.unit   <- beta/sqrt(sum(beta^2))
doPCA.omega <- optim_dopca(beta, Sigma_e, Sigma_XX = Sigma_X)

# TT and Outcome
Z <- (X %*% beta) |> as.numeric()
Y <- (2*Z + 5*C[,1] + rnorm(n)) |> as.numeric()

# Solution
est      <- estimate_DR_curve(Y, Z, C, Z)
solution <- doPCA(Y, X, C, maxit = 1000)

# Visualize objective -----------------------------------------------------
f <- function(omega, Y, X, C){
  print(omega)
  Z <- as.numeric( X %*% omega )
  DR <- estimate_DR_curve(Y, Z, C, Z)$DR_curve
  var(DR)
}

# Theta in (-pi, pi). Evaluate objective for 2D treatment.
theta      <- seq(-pi, pi, length.out = 500)
omega_grid <- cbind(cos(theta), sin(theta))
vals       <- apply(omega_grid, 1, f, Y, X, C)

# For a given problem, plot the objective function as value of theta (which serves as approximation to omega)
png("results/likelihood_function.png", width = 800, height = 600)
plot(theta, vals, type = "l", 
     xlab = expression(theta), 
     ylab = "Objective",
     main = "doPCA Objective over Unit Circle")

# doPCA
abline(v = get_theta(omega), col = "red")
abline(v = get_theta(solution$omega))
abline(v = get_theta(doPCA.omega), col = "blue")

dev.off()