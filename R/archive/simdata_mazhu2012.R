library(MASS)                # mvrnorm
logit_inv <- function(z) 1/(1+exp(-z))
make_sigma <- function(p) toeplitz(0.5^(0:(p-1)))   # Σ_ij = 0.5^{|i–j|}

# beta helpers (pad with zeros for p = 12)
beta1_vec <- function(p) {c(rep(1,6), rep(0, max(0,p-6))) / sqrt(6)}
beta2_vec <- function(p) {c(1,-1,1,-1,1,-1, rep(0, max(0,p-6))) / sqrt(6)}

# predictor generators -------------------------------------------------
gen_X_case1 <- function(n, p){
  if(!p %in% c(6,12))
    stop("Case 1 defined only for p = 6 or 12")
  X <- matrix(NA_real_, n, p)
  if(p==6){
    base <- mvrnorm(n, mu=rep(0,2), Sigma=make_sigma(2))
    X[,1:2] <- base
  } else {                        # p = 12
    base <- mvrnorm(n, mu=rep(0,8), Sigma=make_sigma(8))
    X[,c(1,2,7:12)] <- base
  }
  X
}
gen_X_case2 <- function(n, p){
  mvrnorm(n, mu=rep(0,p), Sigma=make_sigma(p))
}
gen_X <- function(n, p, case){
  if(case==1) gen_X_case1(n,p) else gen_X_case2(n,p)
}

# ----------------------------------------------------------------------
# master builder used by all four models
# ----------------------------------------------------------------------
build_data <- function(n, p, case, y_formula){
  X <- gen_X(n, p, case)
  if(case==1){
    e1 <- rnorm(n); e2 <- rnorm(n)
    X[,3] <- abs(X[,1]+X[,2])         + abs(X[,1]) * e1
    X[,4] <- abs(X[,1]+X[,2])^2       + abs(X[,2]) * e2
    X[,5] <- rbinom(n,1, logit_inv(X[,2]))
    X[,6] <- rbinom(n,1, pnorm(X[,2]))
  }
  beta1 <- beta1_vec(p)
  beta2 <- beta2_vec(p)
  eta1  <- as.numeric(X %*% beta1)
  eta2  <- as.numeric(X %*% beta2)
  Y     <- y_formula(eta1, eta2)
  list(X = X, Y = Y,
       beta1 = beta1, beta2 = beta2,
       eta1 = eta1, eta2 = eta2)
}

# ----------------------------------------------------------------------
# Model-1   𝙔 = η₁ / (0.5 + (η₂+1.5)²) + 0.5 ε
# ----------------------------------------------------------------------
sim_model1 <- function(n, p = 6, case = 1, seed = NULL){
  if(!is.null(seed)) set.seed(seed)
  build_data(
    n, p, case,
    y_formula = function(eta1, eta2){
      eta1 / (0.5 + (eta2 + 1.5)^2) + 0.5 * rnorm(length(eta1))
    })
}

# ----------------------------------------------------------------------
# Model-2   𝙔 = η₁² + 2|η₂| + 0.1|η₂| ε
# ----------------------------------------------------------------------
sim_model2 <- function(n, p = 6, case = 1, seed = NULL){
  if(!is.null(seed)) set.seed(seed)
  build_data(
    n, p, case,
    y_formula = function(eta1, eta2){
      eta1^2 + 2*abs(eta2) + 0.1*abs(eta2)*rnorm(length(eta1))
    })
}

# ----------------------------------------------------------------------
# Model-3   𝙔 = exp(η₁) + 2(η₂+1)² + |η₁| ε
# ----------------------------------------------------------------------
sim_model3 <- function(n, p = 6, case = 1, seed = NULL){
  if(!is.null(seed)) set.seed(seed)
  build_data(
    n, p, case,
    y_formula = function(eta1, eta2){
      exp(eta1) + 2*(eta2 + 1)^2 + abs(eta1)*rnorm(length(eta1))
    })
}

# ----------------------------------------------------------------------
# Model-4   𝙔 = η₁² + η₂² + 0.5 ε
# ----------------------------------------------------------------------
sim_model4 <- function(n, p = 6, case = 1, seed = NULL){
  if(!is.null(seed)) set.seed(seed)
  build_data(
    n, p, case,
    y_formula = function(eta1, eta2){
      eta1^2 + eta2^2 + 0.5*rnorm(length(eta1))
    })
}
