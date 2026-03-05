library(testthat)
library(SuperLearner)
library(mvtnorm)
library(MAVE)

source("R/csMAVE.R")
source("R/compute_new_response_and_exposure.R")
source("R/nuisance_outcome_regression.R")
source("R/nuisance_outcome_C.R")
source("R/nuisance_gps.R")
source("R/estimate_ERS.R")
source("R/estimate_pseudo_outcomes.R")
source("R/estimate_residualized_pair.R")

# ---------------------------------------------------------
# Define Lightweight Fitters for Fast Testing
# ---------------------------------------------------------
fast_lm_fitter <- function(target_var, predictors_df, ...) {
  df <- cbind(target_var = target_var, predictors_df)
  lm(target_var ~ ., data = df)
}

# ---------------------------------------------------------
# Test Suite
# ---------------------------------------------------------

test_that("csMAVE correctly executes the end-to-end SDR pipeline", {
  
  # ---------------------------------------------------------
  # 1. Setup Data with a Known Single-Index Structure
  # ---------------------------------------------------------
  set.seed(123)
  n <- 100
  p <- 3
  q <- 2
  
  # Confounders
  C <- matrix(rnorm(n * q), n, q)
  colnames(C) <- c("C1", "C2")
  
  # Treatments (Confounded by C)
  X <- C %*% matrix(c(0.5, -0.2, 0.1, -0.5, 0.4, 0.8), nrow = 2, ncol = p) + 
    matrix(rnorm(n * p), n, p)
  colnames(X) <- c("X1", "X2", "X3")
  
  # Outcome depends on X ONLY through the index (X1 + X2) and depends on C
  index <- X[, 1] + X[, 2] 
  Y <- as.numeric(C[, 1] - C[, 2] + index + rnorm(n, sd = 0.5))
  
  # =======================================================
  # Test 1: Method 'RA' (Standard exposure, new response)
  # =======================================================
  suppressWarnings({
    res_ra <- csMAVE(Y = Y, X = X, C = C, 
                     method = "RA",
                     args_compute_new_response = list(
                       L = 2, # Small L for speed
                       outcome_fitter = fast_lm_fitter
                     ))
  })
  
  # A. Check Top-Level Structure
  expect_type(res_ra, "list")
  expect_named(res_ra, c("mave_fit", "mave_dim_obj", "d_hat", "new_data", "metadata"))
  
  # B. Verify MAVE Objects
  expect_s3_class(res_ra$mave_fit, "mave")
  expect_s3_class(res_ra$mave_dim_obj, "mave.dim")
  expect_true(is.numeric(res_ra$d_hat) || is.na(res_ra$d_hat))
  
  # C. Verify Data Provenance
  expect_length(res_ra$new_data$new_Y, n)
  expect_equal(dim(res_ra$new_data$new_X), c(n, p))
  
  # Because RA does not transform X, the new_X should equal the original X
  expect_equal(as.numeric(res_ra$new_data$new_X), as.numeric(X))
  
  # =======================================================
  # Test 2: Method 'RP' (Modifies both exposure and response)
  # =======================================================
  suppressWarnings({
    res_rp <- csMAVE(Y = Y, X = X, C = C, 
                     method = "RP",
                     args_compute_new_response = list(
                       L = 2, 
                       C_fitter = fast_lm_fitter
                     ))
  })
  
  expect_s3_class(res_rp$mave_fit, "mave")
  
  # Because RP residualizes X, new_X MUST NOT equal the original X
  expect_false(isTRUE(all.equal(as.numeric(res_rp$new_data$new_X), as.numeric(X))))
  
  # =======================================================
  # Test 3: Verify modifyList Argument Overriding
  # =======================================================
  # We test passing a custom MAVE method string to ensure utils::modifyList 
  # successfully overwrites the default "meanOPG" without crashing.
  
  suppressWarnings({
    res_custom <- csMAVE(Y = Y, X = X, C = C, 
                         method = "RA",
                         args_compute_new_response = list(
                           L = 2, 
                           outcome_fitter = fast_lm_fitter
                         ),
                         args_MAVE = list(method = "meanMAVE"))
  })
  
  # Check that the metadata logged the successful override
  expect_equal(res_custom$metadata$mave_method, "meanMAVE")
  
  # Check that the actual MAVE object successfully used the overridden method
  # (MAVE stores its method call internally)
  expect_true(grepl("MEANMAVE", res_custom$mave_fit$method))
})

test_that("csMAVE correctly runs SuperLearner and recovers the true central subspace", {
  
  # ---------------------------------------------------------
  # 1. Setup Data with a Known Single-Index Structure
  # ---------------------------------------------------------
  set.seed(42)
  n <- 300  # Slightly larger to ensure SuperLearner and MAVE have enough data
  p <- 3
  q <- 2
  
  # Confounders
  C <- matrix(rnorm(n * q), n, q)
  colnames(C) <- c("C1", "C2")
  
  # Treatments (Confounded by C)
  beta_C_to_X <- matrix(c(0.5, -0.2, 0.1, -0.5, 0.4, 0.8), nrow = 2, ncol = p)
  X <- C %*% beta_C_to_X + rmvnorm(n, sigma = diag(p))
  colnames(X) <- c("X1", "X2", "X3")
  
  # TRUE CENTRAL SUBSPACE: Only X1 and X2 matter, equally weighted.
  # The true direction vector must be normalized to length 1 for comparison.
  true_beta_raw <- c(1, 1, 0)
  true_beta <- true_beta_raw / sqrt(sum(true_beta_raw^2))
  
  # Outcome depends on X ONLY through the index (X1 + X2) and depends on C
  index <- X %*% true_beta_raw
  Y <- as.numeric(C[, 1] - C[, 2] + 2 * index + rnorm(n, sd = 0.5))
  
  # ---------------------------------------------------------
  # Common SuperLearner Arguments for Fast Testing
  # ---------------------------------------------------------
  # We use SL.glm to keep the test fast, but it fully exercises the SL wrappers
  test_sl_args <- list(SL.lib = c("SL.glm", "SL.gam"))
  
  # =======================================================
  # Test 1: Method 'RA' (SuperLearner Outcome)
  # =======================================================
  suppressWarnings({
    res_ra <- csMAVE(Y = Y, X = X, C = C, 
                     method = "RA",
                     args_compute_new_response = list(
                       L = 2, 
                       outcome_fitter = SL_outcome_fitter
                     ),
                     args_outcome = test_sl_args)
  })
  
  # Extract the estimated 1D direction from MAVE
  # (MAVE returns a list of matrices in $dir; [[1]] is the 1D subspace)
  est_beta_ra <- res_ra$mave_fit$dir[[1]]
  
  # Calculate Subspace Alignment (Absolute Dot Product)
  alignment_ra <- abs(sum(true_beta * est_beta_ra))
  
  # Expect alignment to be very close to 1 (perfect alignment)
  expect_gt(alignment_ra, 0.90)
  
  # =======================================================
  # Test 2: Method 'DR' (SuperLearner Outcome + SuperLearner GPS)
  # =======================================================
  suppressWarnings({
    res_dr <- csMAVE(Y = Y, X = X, C = C, 
                     method = "DR",
                     args_compute_new_response = list(
                       L = 2, 
                       outcome_fitter = SL_outcome_fitter,
                       gps_fitter = mvn_fitter
                     ),
                     args_outcome = test_sl_args,
                     # Force mvn_fitter to use SuperLearner internally!
                     args_gps = list(method_gps = "SuperLearner", SL.lib = "SL.glm"))
  })
  
  est_beta_dr <- res_dr$mave_fit$dir[[1]]
  alignment_dr <- abs(sum(true_beta * est_beta_dr))
  
  # DR should also successfully recover the true causal subspace
  expect_gt(alignment_dr, 0.90)
  
  # =======================================================
  # Test 3: Method 'RP' (SuperLearner Nuisance C)
  # =======================================================
  suppressWarnings({
    res_rp <- csMAVE(Y = Y, X = X, C = C, 
                     method = "RP",
                     args_compute_new_response = list(
                       L = 2, 
                       C_fitter = SL_nuisance_fitter
                     ),
                     args_C = test_sl_args)
  })
  
  est_beta_rp <- res_rp$mave_fit$dir[[1]]
  alignment_rp <- abs(sum(true_beta * est_beta_rp))
  
  # RP residualizes the exposure entirely. Since the true DGP follows additive 
  # confounding, FWL theorem holds, and it should perfectly recover the subspace.
  expect_gt(alignment_rp, 0.90)
  
  # =======================================================
  # Verify Structural Integrity
  # =======================================================
  expect_s3_class(res_dr$mave_fit, "mave")
  
  # Verify new_data outputs
  expect_length(res_dr$new_data$new_Y, n)
  expect_equal(dim(res_dr$new_data$new_X), c(n, p))
  
  # Verify modifyList safely transported the SL instructions
  expect_equal(res_dr$metadata$cre_pipeline$method, "DR")
})

test_that("csMAVE correctly runs SuperLearner and recovers the true central subspace for all methods", {
  
  # ---------------------------------------------------------
  # 1. Setup Data with a Known Single-Index Structure
  # ---------------------------------------------------------
  set.seed(42)
  n <- 300  # Slightly larger to ensure SuperLearner and MAVE have enough data
  p <- 3
  q <- 2
  
  # Confounders
  C <- matrix(rnorm(n * q), n, q)
  colnames(C) <- c("C1", "C2")
  
  # Treatments (Confounded by C)
  beta_C_to_X <- matrix(c(0.5, -0.2, 0.1, -0.5, 0.4, 0.8), nrow = 2, ncol = p)
  X <- C %*% beta_C_to_X + rmvnorm(n, sigma = diag(p))
  colnames(X) <- c("X1", "X2", "X3")
  
  # TRUE CENTRAL SUBSPACE: Only X1 and X2 matter, equally weighted.
  # The true direction vector must be normalized to length 1 for comparison.
  true_beta_raw <- c(1, 1, 0)
  true_beta <- true_beta_raw / sqrt(sum(true_beta_raw^2))
  
  # Outcome depends on X ONLY through the index (X1 + X2) and depends on C
  index <- X %*% true_beta_raw
  Y <- as.numeric(C[, 1] - C[, 2] + 2 * index + rnorm(n, sd = 0.5))
  
  # ---------------------------------------------------------
  # Common SuperLearner Arguments for Fast Testing
  # ---------------------------------------------------------
  # We use SL.glm to keep the test fast, but it fully exercises the SL wrappers
  test_sl_args <- list(SL.lib = c("SL.glm", "SL.gam"))
  
  # =======================================================
  # Test 1: Method 'RA' (SuperLearner Outcome)
  # =======================================================
  suppressWarnings({
    res_ra <- csMAVE(Y = Y, X = X, C = C, 
                     method = "RA",
                     args_compute_new_response = list(
                       L = 2, 
                       outcome_fitter = SL_outcome_fitter
                     ),
                     args_outcome = test_sl_args)
  })
  
  est_beta_ra <- res_ra$mave_fit$dir[[1]]
  alignment_ra <- abs(sum(true_beta * est_beta_ra))
  
  expect_gt(alignment_ra, 0.90)
  expect_s3_class(res_ra$mave_fit, "mave")
  
  # =======================================================
  # Test 2: Method 'DR' (SuperLearner Outcome + SuperLearner GPS)
  # =======================================================
  suppressWarnings({
    res_dr <- csMAVE(Y = Y, X = X, C = C, 
                     method = "DR",
                     args_compute_new_response = list(
                       L = 2, 
                       outcome_fitter = SL_outcome_fitter,
                       gps_fitter = mvn_fitter
                     ),
                     args_outcome = test_sl_args,
                     args_gps = list(method_gps = "SuperLearner", SL.lib = "SL.glm"))
  })
  
  est_beta_dr <- res_dr$mave_fit$dir[[1]]
  alignment_dr <- abs(sum(true_beta * est_beta_dr))
  
  expect_gt(alignment_dr, 0.90)
  expect_s3_class(res_dr$mave_fit, "mave")
  
  # =======================================================
  # Test 3: Method 'PO' (Pseudo-Outcomes with SuperLearner)
  # =======================================================
  suppressWarnings({
    res_po <- csMAVE(Y = Y, X = X, C = C, 
                     method = "PO",
                     args_compute_new_response = list(
                       L = 2, 
                       outcome_fitter = SL_outcome_fitter,
                       gps_fitter = mvn_fitter
                     ),
                     args_outcome = test_sl_args,
                     args_gps = list(method_gps = "SuperLearner", SL.lib = "SL.glm"))
  })
  
  est_beta_po <- res_po$mave_fit$dir[[1]]
  alignment_po <- abs(sum(true_beta * est_beta_po))
  
  expect_gt(alignment_po, 0.90)
  expect_s3_class(res_po$mave_fit, "mave")
  expect_equal(res_po$metadata$cre_pipeline$method, "PO")
  
  # =======================================================
  # Test 4: Method 'RP' (SuperLearner Nuisance C)
  # =======================================================
  suppressWarnings({
    res_rp <- csMAVE(Y = Y, X = X, C = C, 
                     method = "RP",
                     args_compute_new_response = list(
                       L = 2, 
                       C_fitter = SL_nuisance_fitter
                     ),
                     args_C = test_sl_args)
  })
  
  est_beta_rp <- res_rp$mave_fit$dir[[1]]
  alignment_rp <- abs(sum(true_beta * est_beta_rp))
  
  # RP residualizes the exposure entirely. Since the true DGP follows additive 
  # confounding, FWL theorem holds, and it should perfectly recover the subspace.
  expect_gt(alignment_rp, 0.90)
  expect_s3_class(res_rp$mave_fit, "mave")
  
  # =======================================================
  # Verify Structural Data Integrity (Using PO as the proxy)
  # =======================================================
  
  # Verify new_data outputs for the PO method
  expect_length(res_po$new_data$new_Y, n)
  expect_equal(dim(res_po$new_data$new_X), c(n, p))
  
  # Verify that for PO, the exposure matrix X is NOT modified
  expect_equal(as.numeric(res_po$new_data$new_X), as.numeric(X))
  
  # Verify that for RP, the exposure matrix X IS modified
  expect_false(isTRUE(all.equal(as.numeric(res_rp$new_data$new_X), as.numeric(X))))
})

test_that("csMAVE correctly handles single and multiple methods", {
  
  # 1. Setup Data
  set.seed(42)
  n <- 200 
  p <- 3
  q <- 2
  
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)
  
  C <- matrix(rnorm(n * q), n, q)
  colnames(C) <- paste0("C", 1:q)
  
  Y <- as.numeric(C %*% c(0.5, -0.5) + X %*% c(1, -1, 0) + rnorm(n))
  
  # =======================================================
  # Test 1: Single Method (Backwards Compatibility)
  # =======================================================
  suppressWarnings({
    res_single <- csMAVE(Y = Y, X = X, C = C, 
                         method = "DR",
                         args_compute_new_response = list(
                           L = 2,
                           outcome_fitter = fast_lm_fitter,
                           gps_fitter = mvn_fitter
                         ))
  })
  
  expect_type(res_single, "list")
  expect_s3_class(res_single$mave_fit, "mave")
  expect_true("d_hat" %in% names(res_single))
  expect_equal(res_single$metadata$causal_method, "DR")
  
  # =======================================================
  # Test 2: Multiple Methods Simultaneously
  # =======================================================
  suppressWarnings({
    res_multi <- csMAVE(Y = Y, X = X, C = C, 
                        method = c("DR", "PO", "RP"),
                        args_compute_new_response = list(
                          L = 2,
                          outcome_fitter = fast_lm_fitter,
                          gps_fitter = mvn_fitter,
                          C_fitter = fast_lm_fitter
                        ))
  })
  
  expect_type(res_multi, "list")
  # The top level should now be a named list of the methods
  expect_named(res_multi, c("DR", "PO", "RP"))
  
  # Each sub-element should be a fully formed MAVE object
  expect_s3_class(res_multi$DR$mave_fit, "mave")
  expect_s3_class(res_multi$PO$mave_fit, "mave")
  expect_s3_class(res_multi$RP$mave_fit, "mave")
  
  expect_equal(res_multi$DR$metadata$causal_method, "DR")
  expect_equal(res_multi$PO$metadata$causal_method, "PO")
})
