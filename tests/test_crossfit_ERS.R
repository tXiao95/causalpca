library(testthat)
library(data.table)

# Source your functions
source("R/crossfit_ERS.R")
source("R/nuisance_outcome_regression.R")
source("R/nuisance_gps.R")

# ---------------------------------------------------------
# Define Lightweight Fitters for Fast Testing
# ---------------------------------------------------------
# NEW TEST FITTER (Put this in test_crossfit_ERS.R)
fast_lm_fitter <- function(Y, XC_df, ...) {
  df <- cbind(Y = Y, XC_df)
  lm(Y ~ ., data = df)
}

# ---------------------------------------------------------
# Test Suite
# ---------------------------------------------------------
test_that("crossfit_ERS correctly evaluates globally and handles L=1 vs L>1", {
  
  # 1. Setup Data
  set.seed(42)
  n <- 300 
  p <- 2
  q <- 2
  
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- c("X1", "X2")
  
  C <- matrix(rnorm(n * q), n, q)
  colnames(C) <- c("C1", "C2")
  
  Y <- as.numeric(C %*% c(0.5, -0.5) + X %*% c(1, -1) + rnorm(n))
  
  # Define a small evaluation grid
  x_eval <- rbind(c(0, 0), c(1, -1), c(-1, 1))
  colnames(x_eval) <- c("X1", "X2")
  
  # =======================================================
  # Test 1: L=1 (Behavior identical to old estimate_ERS)
  # =======================================================
  suppressWarnings({
    res_L1 <- crossfit_ERS(Y = Y, X = X, C = C, x_eval = x_eval,
                           estimator = "DR",
                           L = 1,
                           outcome_fitter = SL_outcome_fitter,
                           gps_fitter = mvn_fitter,
                           args_gps = list(method_gps = "linear"),
                           args_outcome = list(SL.lib = c("SL.glm", "SL.earth", "SL.xgboost")), 
                                               optimize_bw = TRUE)
  })
  
  expect_type(res_L1, "list")
  expect_named(res_L1, c("results", "metadata"))
  expect_equal(res_L1$metadata$L_folds, 1)
  expect_equal(nrow(res_L1$results), 3) # Should match the 3 rows in x_eval
  
  # Ensure the point estimates are numeric and non-NA
  expect_true(all(!is.na(res_L1$results$estimate)))
  expect_true(all(!is.na(res_L1$results$se)))
  
  # =======================================================
  # Test 2: L=5 (Cross-Fitting with Global Integration)
  # =======================================================
  suppressWarnings({
    res_L5 <- crossfit_ERS(Y = Y, X = X, C = C, x_eval = x_eval,
                           estimator = "DR",
                           L = 5,
                           outcome_fitter = fast_lm_fitter,
                           gps_fitter = mvn_fitter,
                           args_gps = list(method_gps = "linear"), optimize_bw = TRUE)
  })
  
  expect_type(res_L5, "list")
  expect_equal(res_L5$metadata$L_folds, 5)
  expect_equal(nrow(res_L5$results), 3)
  
  # Cross-fitted estimates should be reasonably close to the L=1 estimates
  # (They won't be perfectly identical due to fold splitting, but should be highly correlated)
  expect_true(cor(res_L1$results$estimate, res_L5$results$estimate) > 0.9)
  
  # =======================================================
  # Test 3: Regression Adjustment (RA) Fallback
  # =======================================================
  suppressWarnings({
    res_ra <- crossfit_ERS(Y = Y, X = X, C = C, x_eval = x_eval,
                           estimator = "RA",
                           L = 3,
                           outcome_fitter = fast_lm_fitter)
  })
  
  expect_equal(res_ra$metadata$estimator, "RA")
  expect_true(all(is.na(res_ra$results$se))) # RA doesn't compute standard errors in this implementation
  expect_true(all(!is.na(res_ra$results$estimate)))
  
  # =======================================================
  # Test 4: Bandwidth Optimization (optimize_bw = TRUE)
  # =======================================================
  suppressWarnings({
    res_opt <- crossfit_ERS(Y = Y, X = X, C = C, x_eval = x_eval,
                            estimator = "DR",
                            L = 2,
                            outcome_fitter = fast_lm_fitter,
                            gps_fitter = mvn_fitter,
                            optimize_bw = TRUE,
                            args_gps = list(method_gps = "linear"))
  })
  
  expect_true(res_opt$metadata$optimized_bw)
  # Check that the dynamic bandwidths were assigned and are numeric
  expect_true(all(!is.na(res_opt$results$h_used_1)))
  expect_type(res_opt$results$h_used_1, "double")
  
  # Because undersmoothing is applied in optimize_bw, the final bandwidths 
  # should generally be different from the default h_pilot
  expect_false(isTRUE(all.equal(res_L5$results$h_used_1, res_opt$results$h_used_1)))
})
