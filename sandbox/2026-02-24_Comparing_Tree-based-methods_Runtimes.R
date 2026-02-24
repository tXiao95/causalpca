library(ggplot2)
library(here)
library(microbenchmark)
library(ranger)
library(earth)
library(xgboost)
library(lightgbm)
library(dplyr)
library(tidyr)

# Helper to generate synthetic data with 20 covariates
generate_data <- function(n, p = 20) {
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)
  # Non-linear outcome with interactions
  Y <- 2*X[,1] + 1.5*X[,2]^2 + X[,3]*X[,4] + rnorm(n)
  return(list(X = X, Y = Y))
}

# Define the sample sizes to test (logarithmic scale)
n_sizes <- c(1000, 5000, 10000, 25000, 50000, 100000)
#n_sizes <- c(1000, 5000)
results_list <- list()

for (n in n_sizes) {
  message("Benchmarking n = ", n)
  data <- generate_data(n)
  
  # Standardize XGBoost/LightGBM inputs
  dtrain <- xgb.DMatrix(data = data$X, label = data$Y)
  lgb_train <- lgb.Dataset(data = data$X, label = data$Y)
  
  bm <- microbenchmark(
    Ranger = {
      ranger(y = data$Y, x = data$X, num.trees = 100, num.threads = 1)
    },
    Earth = {
      earth(x = data$X, y = data$Y)
    },
    XGBoost = {
      xgb.train(params = list(objective = "reg:squarederror", nthread = 1), 
                data = dtrain, nrounds = 100)
    },
    LightGBM = {
      lgb.train(params = list(objective = "regression", num_threads = 1, 
                              learning_rate = 0.1, verbose = -1), 
                data = lgb_train, nrounds = 100)
    },
    times = 3 # Low reps because large n is slow
  )
  
  df_n <- as.data.frame(bm)
  df_n$n <- n
  results_list[[as.character(n)]] <- df_n
}

# Combine and process results
res_df <- bind_rows(results_list) %>%
  mutate(time_sec = time / 1e9) %>% # Convert nanoseconds to seconds
  group_by(n, expr) %>%
  summarise(mean_time = mean(time_sec),
            sd_time = sd(time_sec), .groups = 'drop')

# 3. Plotting the results
plt <- ggplot(res_df, aes(x = n, y = mean_time, color = expr)) +
  geom_line(linewidth = 1) +
  geom_point() +
  scale_x_continuous(labels = scales::comma) +
  labs(title = "Algorithm Runtime Comparison (Standard Regression)",
       subtitle = "Complexity scaling with n (20 covariates)",
       x = "Sample Size (n)",
       y = "Mean Runtime (Seconds)",
       color = "Algorithm") +
  theme_minimal()

resultspath <- here("outputs", "experiments", "2026-02-24_Comparing_Tree-based-methods_runtime.pdf")
resultspathtab <- here("outputs", "experiments", "2026-02-24_Comparing_Tree-based-methods_runtime.csv")


readr::write_csv(res_df, resultspathtab)
ggsave(resultspath)
