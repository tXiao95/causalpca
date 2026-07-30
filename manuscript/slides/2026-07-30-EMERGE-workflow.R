# Reproduce the R code and console results shown on workflow slides 13--16
# of 2026-07-30-EMERGE.qmd.

library(csdr)
data(csdr_example2)

Y <- csdr_example2$Y
A <- csdr_example2$A
C <- csdr_example2$C

head(Y, 5)
head(A[, c(1:3, 18:20)], 5)
head(C, 5)

sl_lib <- c(
  "SL.glm", "SL.earth", "SL.ranger"
)

learners <- csdr_learners(
  sl_library = sl_lib
)

set.seed(42)
fit <- csdr(
  Y = Y, A = A, C = C,
  L = 5, seed = 42,
  learners = learners,
  max_dim = 5,
  keep_nuisance = TRUE
)

summary(fit)

B_hat <- coef(fit)
Z_hat <- scores(fit)

mave <- mave_fits(fit)
dim_selection <- mave$dimension_selection

B_hat[c(1, 2, 20), , drop = FALSE]
sis(fit)
head(Z_hat, 2)
dim_selection

z_grid <- matrix(
  rep(apply(Z_hat, 2, median), each = 50),
  nrow = 50,
  dimnames = list(NULL, colnames(Z_hat))
)
z_grid[, 1] <- seq(
  quantile(Z_hat[, 1], 0.05),
  quantile(Z_hat[, 1], 0.95),
  length.out = 50
)

ers <- estimate_ers(
  Y = Y, A = Z_hat, C = C,
  a_eval = z_grid, estimator = "DR",
  L = 5, seed = 42,
  verbose = FALSE
)

head(ers$results)
