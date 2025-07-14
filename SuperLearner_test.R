
set.seed(1)
source("R/doPCA.R")
source("R/estimate_DR_curve.R")

# Settings
n <- 200; p <- 2; q <- 2
A <- matrix(rnorm(p * q, 1, 1), nrow = p, ncol = q)
C <- matrix(rnorm(n*q, mean = 1), n, q)
X <- C %*% t(A) + rnorm(n*p)

# Quantities
Sigma_C
Sigma_XX <- Adiag(1, nrow = p)

omega_true <- c(1,0.5); omega_true <- omega_true/sqrt(sum(omega_true^2))

# TT and Outcome
Z <- (X %*% omega_true) |> as.numeric()
Y <- (2*Z + 5*C[,1] + rnorm(n)) |> as.numeric()

# Solution
est      <- estimate_DR_curve(Y, Z, C, Z)
solution <- doPCA(Y, X, C, maxit = 1000)

f <- function(omega, Y, X, C){
  print(omega)
  Z <- as.numeric( X %*% omega )
  DR <- estimate_DR_curve(Y, Z, C, Z)$DR_curve
  var(DR)
}


f(omega_true, Y, X, C)
f(solution$omega, Y, X, C)

theta      <- seq(0, pi, length.out = 250)
omega_grid <- cbind(cos(theta), sin(theta))

vals <- apply(omega_grid, 1, f, Y, X, C)

vals <- vals[1:250]

# For a given problem, plot the objective function as value of theta (which serves as approximation to omega)
png("results/likelihood_function.png", width = 800, height = 600)
plot(theta, vals, type = "l", 
     xlab = expression(theta), 
     ylab = "Objective",
     main = "doPCA Objective over Unit Circle")
abline(v = 0.463648, col = "red")
abline(v = 0.283)

dev.off()


# On fourth iteration: omega is
# (c(0, -7.14180690193068e-17, 0))

# If you want to use equality constraints, then you should use one of these algorithms 
# NLOPT_LD_AUGLAG, NLOPT_LN_AUGLAG, NLOPT_LD_AUGLAG_EQ, NLOPT_LN_AUGLAG_EQ, NLOPT_GN_ISRES, 
# NLOPT_LD_SLSQP, NLOPT_LN_COBYLA

# Experiments -------------------------------------------------------------

# Re-parameterization trick
# Equality constraint + derivative-free 
# Calculate first and second derivatives numerically and use Newton-Raphson
# Other solutions? 

# Need to visualize my objective function


# Linear model
# Nonlinear
# Interactions