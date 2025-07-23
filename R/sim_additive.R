library(MASS)

#' @param n Integer — sample size.
#' @param p Integer — dimension of \eqn{X}.  Default \code{6}.
#' @param q Integer — dimension of \eqn{C}.  Default \code{4}.
#' @param beta Numeric length-\code{p} vector.  If \code{NULL} (default)
#'   the vector \code{c(1, 1, 1, 0, 0,...)} is used.
#' @param gamma Numeric length-\code{q} vector.  If \code{NULL} (default)
#'   a vector of ones is created and scaled to unit Euclidean length.
#' @param A \code{p × q} loading matrix linking \eqn{C} to \eqn{X}.  If
#'   \code{NULL}, a tiled identity matrix is built so that \emph{every}
#'   row of \eqn{A} contains at least one non-zero entry (even when
#'   \code{p > q}).
#' @param f Unary function implementing \(f\).  Default is the identity.
#' @param g Unary function implementing \(g\).  Default is the identity.
#' @param e.dist Character string: distribution of the measurement error
#'   \eqn{e}.  One of \code{"normal"}, \code{"uniform"} or
#'   \code{"laplace"}.  Default \code{"normal"}.
#' @param e.scale Positive numeric scale parameter for \eqn{e}
#'   (sd for Normal, half-range for Uniform, Laplace scale parameter).
#'   Default \code{1}.
#' @param eps.sd Standard deviation of the additive noise
#'   \eqn{\varepsilon}.  Default \code{1}.
#' @param seed Optional integer seed for reproducibility.
sim_additive <- function(n,
                         p = 6,
                         q = 4,
                         beta  = NULL,   # default (1,1,1,0,0,…)
                         gamma = NULL,   # default all-ones/√q
                         lambda = NULL,  # default 0
                         A     = NULL,   # p×q loading matrix from C to X
                         f = function(z) z,
                         g = function(z) z,
                         e.dist  = "normal",   # "normal", "uniform", "laplace"
                         e.scale = 2,
                         eps.sd  = 1,
                         confounding_type = "additive",
                         seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  ## default beta / gamma --------------------------------------------------
  if (is.null(beta)) {
    beta <- numeric(p); beta[1:min(3, p)] <- 1      # (1,1,0,…)
  }
  if (is.null(gamma)) {
    gamma <- rep(1, q) / sqrt(q)
  }
  if (is.null(lambda)) {
    lambda <- 0
  }
  stopifnot(length(beta) == p, length(gamma) == q)
  
  ## default A -------------------------------------------------------------
  if (is.null(A)) {
    ## Idea: repeat the q×q identity downward until we have p rows,
    ## then keep the first p rows.  Each X_j is linked to exactly
    ## one C_k, and every C_k affects roughly ⌈p/q⌉ components of X.
    reps <- ceiling(p / q)
    A <- diag(q)[ rep(seq_len(q), reps), , drop = FALSE ][1:p, , drop = FALSE]
  }
  stopifnot(nrow(A) == p, ncol(A) == q)
  
  ## draw e ----------------------------------------------------------------
  draw_e <- function() {
    dist <- match.arg(e.dist, c("normal", "uniform", "laplace"))
    nd   <- n * p
    vec <- switch(dist,
                  normal  = rnorm(nd, sd = e.scale),
                  uniform = runif(nd, -e.scale,  e.scale),
                  laplace = {u <- runif(nd, -0.5, 0.5)
                  e.scale * sign(u) * log1p(-2 * abs(u))})
    matrix(vec, n, p)
  }
  
  ## generate data ---------------------------------------------------------
  C <- matrix(rnorm(n * q), n, q) + 1 # N(1_p, I_p)
  
  if(confounding_type == "additive"){
    X <- ( C %*% t(A) ) + draw_e()          # AC + e
  } else if(confounding_type == "multiplicative"){
    X <- ( C %*% t(A) ) * draw_e()          # ACe
  }
  
  eta.x <- as.numeric(X %*% beta)
  eta.c <- as.numeric(C %*% gamma)
  
  Y <- f(eta.x) + g(eta.c) + lambda * (eta.x * eta.c) + rnorm(n, sd = eps.sd)
  
  list(Y = Y, X = X, C = C,
       eta.x = eta.x, eta.c = eta.c,
       beta = beta, gamma = gamma, lambda = lambda, A = A,
       e.dist = e.dist, e.scale = e.scale, eps.sd = eps.sd)
}
