library(here)
library(data.table)
library(ggplot2)

source(here("R", "sims_from_papers.R"))
source(here("R", "estimate_Seff.R"))

run_scaling_experiment <- function(n_reps = 500, p = 12) {
  
  # Define the sequence of sample sizes
  n_seq <- c(500, 1000, 5000, 10000, 50000)
  #n_seq <- c(150000)
  #n_seq <- c(275000)
  
  # 1. Pre-allocate a long-format dataframe for clean aggregation and plotting later
  total_runs <- length(n_seq) * n_reps
  results_df <- data.frame(
    n = integer(total_runs),
    rep = integer(total_runs),
    Method = character(total_runs),
    Frobenius_Dist = numeric(total_runs),
    Runtime_Secs = numeric(total_runs),
    stringsAsFactors = FALSE
  )
  
  cat("Starting scaling experiment for Model IV (p =", p, ")...\n")
  row_idx <- 1
  
  for (n_val in n_seq) {
    cat(sprintf("\n--- Running Sample Size: N = %d ---\n", n_val))
    
    for (i in seq_len(n_reps)) {
      # Dynamic seed for reproducibility
      current_seed <- i * 1000 + n_val
      
      # Generate Model 4 data
      sim_data <- sim_luo_cai_mod4(n = n_val, p = p, seed = current_seed)
      X <- sim_data$X
      Y <- sim_data$Y
      beta_true <- sim_data$beta_true
      d <- sim_data$d
      
      # Calculate Theoretical Bandwidths for EE
      b_val <- n_val^(-1 / (d + 4))
      h_val <- n_val^(-1 / (4 * p))
      
      # ---------------------------------------------------------
      # Method 1: meanMAVE
      # ---------------------------------------------------------
      t0 <- Sys.time()
      suppressWarnings({
        fit_mave <- MAVE::mave(Y ~ X, method = "meanMAVE")
      })
      t1 <- Sys.time()
      beta_mave <- fit_mave$dir[[d]]
      
      results_df[row_idx, ] <- list(
        n_val, i, "meanMAVE", 
        Delta(beta_true, beta_mave, "F"), 
        as.numeric(difftime(t1, t0, units = "secs"))
      )
      row_idx <- row_idx + 1
      
      # ---------------------------------------------------------
      # Method 2: meanOPG
      # ---------------------------------------------------------
      t0 <- Sys.time()
      suppressWarnings({
        fit_opg <- MAVE::mave(Y ~ X, method = "meanOPG")
      })
      t1 <- Sys.time()
      beta_opg <- fit_opg$dir[[d]]
      
      results_df[row_idx, ] <- list(
        n_val, i, "meanOPG", 
        Delta(beta_true, beta_opg, "F"), 
        as.numeric(difftime(t1, t0, units = "secs"))
      )
      row_idx <- row_idx + 1
      
      # ---------------------------------------------------------
      # Method 3: Efficient Score (1 Step)
      # ---------------------------------------------------------
      t0 <- Sys.time()
      beta_ee_1step <- run_efficient_estimator(
        X = X,
        Y = Y,
        beta_init = beta_mave,
        b = b_val,
        h = h_val,
        max_iters = 0, # Force a single Newton-Raphson step
        SL = FALSE
      )
      t1 <- Sys.time()
      
      results_df[row_idx, ] <- list(
        n_val, i, "EE_1Step", 
        Delta(beta_true, beta_ee_1step, "F"), 
        as.numeric(difftime(t1, t0, units = "secs"))
      )
      row_idx <- row_idx + 1
      
      # Print a tiny progress tracker every 10 iterations
      if (i %% 10 == 0) cat(i, " ")
    }
    cat("\n")
  }
  
  cat("\n--- Experiment Complete ---\n")
  return(results_df)
}

# Execute the experiment
# Warning: Monitor cluster memory usage closely for n >= 50,000
scaling_results <- run_scaling_experiment(n_reps = 5, p = 12)
resultspath <- here("outputs", "experiments", "scaling_runtimes_Model4_p12_50K.rds")
resultspathplot <- here("outputs", "experiments", "scaling_runtimes_Model4_p12_50K.pdf")

dt <- data.table(scaling_results)
dt <- dt[, .(dist_F = mean(Frobenius_Dist), 
             dist_F_sd = sd(Frobenius_Dist),
             runtime_min = mean(Runtime_Secs / 60),
             runtime_min_sd = sd(Runtime_Secs / 60)), by = .(n, Method)]

dt <- readRDS(resultspath)

# Plot power law ----------------------------------------------------------

dt <- dt[Method != "EE_1Step"]

# 2. Fit a Quadratic Complexity Model: Runtime = c1*n + c2*n^2
# The "0 +" removes the intercept so n=0 takes exactly 0 minutes
# The interaction ":" ensures meanMAVE and meanOPG get their own unique c1 and c2 coefficients
fit_poly <- lm(runtime_min ~ 0 + Method:n + Method:I(n^2), data = dt)

cat("Quadratic Model Coefficients (c1*n + c2*n^2):\n")
print(coef(fit_poly))

# 3. Create a projection grid of new sample sizes (matching your image's x-axis)
n_project <- seq(500, 275000, length.out = 200)

grid_eval <- expand.grid(
  n = n_project,
  Method = c("meanMAVE", "meanOPG")
)

# 4. Predict the runtimes using the new quadratic model
grid_eval$pred_runtime <- predict(fit_poly, newdata = grid_eval)

# Print some milestones to check the new predictions
milestones <- grid_eval #|> filter(n %in% c(50000, 100000, 250000))
print(milestones)

# 5. Plot the observed data against the projected curves
ggplot() +
  # Projected Curves
  geom_line(data = grid_eval, aes(x = n, y = pred_runtime, color = Method), 
            linewidth = 1, linetype = "dashed") +
  # Observed Data Points
  geom_point(data = dt, aes(x = n, y = runtime_min, color = Method), 
             size = 4) +
  scale_x_continuous(labels = scales::comma_format()) +
  scale_y_continuous(labels = scales::comma_format()) +
  labs(
    title = "Empirical Runtime Projection",
    subtitle = expression("Quadratic fit based on algorithmic complexity: " ~ c[1]*n + c[2]*n^2),
    x = "Sample Size (n)",
    y = "Runtime (Minutes)"
  ) +
  theme_minimal() +
  theme(
    legend.position = c(0.2, 0.8),
    legend.background = element_rect(fill = "white", color = "gray90"),
    text = element_text(size = 12)
  )
# Save the raw data
saveRDS(dt, resultspath)
ggsave(resultspathplot)
