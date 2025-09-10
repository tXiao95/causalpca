
# Ma & Zhu (2013) Annals of Statistics Simulation Examples ----------------

generate_example1 <- function(n, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  beta <- c(1.3, -1.3, 1.0, -0.5, 0.5, -0.5)
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  e1 <- rnorm(n)
  e2 <- rnorm(n)
  X3 <- 0.2*X1 + 0.2*(X2 + 2)^2 + 0.2*e1
  X4 <- 0.1 + 0.1*(X1 + X2) + 0.3*(X1 + 1.5)^2 + 0.2*e2
  invlogit <- function(z) 1/(1+exp(-z))
  X5 <- rbinom(n, 1, invlogit(X1))
  X6 <- rbinom(n, 1, invlogit(X2))
  X <- cbind(X1,X2,X3,X4,X5,X6)
  Y <- as.numeric(X %*% beta) + rnorm(n, sd = 1)
  list(X = X, Y = Y, beta = beta)
}


generate_example2 <- function(n, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  dat <- generate_example1(n)  # reuse X and beta construction
  X    <- dat$X
  beta <- dat$beta
  lin  <- as.numeric(X %*% beta)
  mu   <- sin(2 * lin) + 2 * exp(2 + lin)
  sigma <- sqrt(log(2 + lin^2))
  Y <- rnorm(n, mean = mu, sd = sigma)
  list(X = X, Y = Y, beta = beta)
}

generate_example3 <- function(n, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  beta1 <- c(1, 2/3, 2/3, 0, -1/3, 2/3)
  beta2 <- c(0.8, 0.8, -0.3, 0.3, 0, 0)
  U1 <- sample(c(1, -1), n, replace = TRUE, prob = c(0.5, 0.5))
  U2 <- sample(c(sqrt(3/7), -sqrt(7/3)), n, replace = TRUE, prob = c(0.7, 0.3))
  U3 <- runif(n, -sqrt(3), sqrt(3))
  U4 <- runif(n, -sqrt(3), sqrt(3))
  U5 <- runif(n, -sqrt(3), sqrt(3))
  U6 <- runif(n, -sqrt(3), sqrt(3))
  X1 <- U1 - U2
  X2 <- U2 - U3 - U4
  X3 <- U3 + U4
  X4 <- 2 * U4
  X5 <- U5 + 0.5 * U6
  X6 <- U6
  X <- cbind(X1,X2,X3,X4,X5,X6)
  lin1 <- as.numeric(X %*% beta1)
  lin2 <- as.numeric(X %*% beta2)
  mu <- 2 * lin1^2
  sigma <- sqrt(2 * exp(lin2))
  Y <- rnorm(n, mean = mu, sd = sigma)
  list(X = X, Y = Y, beta1 = beta1, beta2 = beta2)
}
