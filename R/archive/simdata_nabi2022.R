simdata_nabi <- function(n, p = 6, case = 1, seed = NULL) {
  if (!p %in% c(6, 12))  stop("p must be 6 or 12 (as in the paper)")
  if (!case %in% c(1, 2)) stop("case must be 1 or 2")
  if (!is.null(seed)) set.seed(seed)
  
  ## ------------------------------------------------------------
  ## helper objects
  ## ------------------------------------------------------------
  ar_cov   <- function(d) toeplitz(0.5 ^ (0:(d - 1)))  # Σ_ij = 0.5^{|i-j|}
  invlogit <- function(z) 1 / (1 + exp(-z))
  
  ## baseline covariates  (same in both cases)
  C     <- matrix(rnorm(n * 4), n, 4)                    # N4(0, I)
  Csum  <- rowSums(C)
  Calt  <- C[, 1] - C[, 2] - C[, 3] + C[, 4]          #  C1 − C2 − C3 + C4
  Csign <- C %*% ((-1)^(1:4))                        #  Σ (-1)^i Ci
  
  ## ------------------------------------------------------------
  ## Treatment A
  ## ------------------------------------------------------------
  A <- matrix(NA_real_, n, p)
  
  if (case == 1) {
    # ---- Case 1: mixed / non-Gaussian -------------------------
    if (p == 6) {
      k <- 2                               # only A1,A2 from MVN
      mu_block <- cbind(Csum, Csign)       # n × 2
      Sigma_block <- ar_cov(k)
      Z <- MASS::mvrnorm(n, mu = rep(0, k), Sigma_block)
      A[, 1:2] <- Z + mu_block
    } else {                               # p == 12
      k <- 8                               # A1, A2, A7:12
      mu_extra <- cbind(C[, 1], C[, 2], C[, 3],
                        -C[, 1] + C[, 2],
                        -C[, 2] + C[, 3],
                        -C[, 3] + C[, 4])  # n × 6
      mu_block <- cbind(Csum, Csign, mu_extra)      # n × 8
      Sigma_block <- ar_cov(k)
      Z <- MASS::mvrnorm(n, mu = rep(0, k), Sigma_block)
      A[, c(1, 2, 7:12)] <- Z + mu_block
    }
    
    ## build the four non-Gaussian coordinates
    e1 <- rnorm(n);  e2 <- rnorm(n)
    A[, 3] <- abs(A[, 1] + A[, 2])              + sqrt(abs(A[, 1])) * e1
    A[, 4] <- abs(A[, 1] + A[, 2])^0.5          + sqrt(abs(A[, 2])) * e2
    A[, 5] <- rbinom(n, 1, invlogit(A[, 2]))
    A[, 6] <- rbinom(n, 1, pnorm(A[, 2]))
  } else {
    # ---- Case 2: multivariate normal --------------------------
    mu <- if (p == 6) {
      cbind(Csum,
            Csign,
            Calt,
            -Calt,
            Csum - 2 * C[, 3],
            Csum - 2 * C[, 1])              # n × 6
    } else {                                # p == 12
      cbind(Csum,
            Csign,
            Calt,
            -Calt,
            Csum - 2 * C[, 3],
            Csum - 2 * C[, 1],
            C[, 1],  C[, 2],  C[, 3],
            -C[, 1], -C[, 2], -C[, 3])      # n × 12
    }
    A <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = ar_cov(p)) + mu
  }
  
  ## ------------------------------------------------------------
  ##   True structural directions & indices
  ## ------------------------------------------------------------
  beta1 <- beta2 <- rep(0, p)
  beta1[1:6] <-  1 / sqrt(6)
  beta2[1:6] <- c(1, -1, 1, -1, 1, -1) / sqrt(6)
  
  eta1 <- as.numeric(A %*% beta1)
  eta2 <- as.numeric((A ^ 2) %*% beta2)           # quadratic index
  
  ## ------------------------------------------------------------
  ##   Outcome Y   (paper’s specification)
  ## ------------------------------------------------------------
  Y <- eta1 + 0*eta2 + Csum + rowSums(A) * Csum + rnorm(n)
  
  ## ------------------------------------------------------------
  ##   return
  ## ------------------------------------------------------------
  list(
    Y     = Y,
    X     = A,
    C     = C,
    beta1 = beta1,
    beta2 = beta2,
    eta1  = eta1,
    eta2  = eta2
  )
}
