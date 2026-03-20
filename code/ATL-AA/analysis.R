library(data.table)
library(ggplot2)
library(here)
library(readxl)

X_path <- here("data/ATL-AA/ATL_AA_PFAS_N=532_two_visits.csv")
C_path <- here("data/ATL-AA/MPTB_AA Cohort Data_N766, clean, DEidentified, data dictionary_12.5.2025.xlsx")

# Data path 
X_dat <- fread(X_path)[Visit == 1, .(subjectid, PFHXS, PFOS, PFOA, PFNA, PFDA, PFUNDA, PFDODA)]
setnames(X_dat, "subjectid", "Subjectid")

C_data <- read_excel(C_path, sheet = 1) |> data.table()
dict   <- read_excel(C_path, sheet = 2) |> data.table()

# Pick covariates
C_dat <- C_data[AllFullTerm == 1, .(Subjectid, 
                    age_enrollment, 
                    Education_4.level,
                #    income_5cat, 
                   # Parity_3cat,
                    FirstPrenatalBMI,
                  #  MarriedCohab_Not,
                    TobaccoUse_MRorSR,
                    #AlcoholUse_MRorSR,
                    MarijuanaUse_MRorSR,
                    Sex,
                    #birthga,
                    birth_weight)]

C_dat[, birth_weight := ifelse(birth_weight == "NA", NA, as.numeric(birth_weight) )]

dat <- merge(X_dat, C_dat, by = "Subjectid")


# Apply log2 transformation to all PFAS concentrations to reduce influenc eof outliers
dat[, `:=`(PFHXS = log2(PFHXS), 
           PFOS = log2(PFOS),
           PFOA = log2(PFOA),
           PFNA = log2(PFNA))]
           #PFDA = log2(PFDA),
           #PFUNDA = log2(PFUNDA),
           #PFDODA = log2(PFDODA))]

# 
dat[, `:=`(edu = factor(Education_4.level, levels = 1:4, labels = c("Less than HS", 
                                                                    "HS or GED", 
                                                                    "Some college or tech school",
                                                                    "4-yr college or more")),
           bmi = FirstPrenatalBMI,
           age = age_enrollment,
           tobacco = factor(TobaccoUse_MRorSR, levels = 0:1, labels = c("No", "Yes")),
           marijuana = factor(MarijuanaUse_MRorSR, levels = 0:1, labels = c("No", "Yes")),
           sex = as.factor(Sex))]

dat <- dat[, .(PFOS, PFOA, PFNA, PFHXS, 
               #PFDA, PFUNDA, PFDODA, 
               age, bmi, edu, tobacco, marijuana, sex, birth_weight)]

# Setup & Libraries -------------------------------------------------------
library(SuperLearner)
library(MAVE)
library(torch) # Required for the Neural Network

# Ensure you have all wrappers sourced
source(here("R/nuisance_outcome_regression.R")) 
source(here("R/nuisance_gps.R")) # Contains mvn_fitter
source(here("R/crossfit_ERS.R")) # Contains crossfit_ERS

X_vars <- c("PFOS", "PFOA", "PFNA", "PFHXS")
C_vars <- c("age", "bmi", "edu", "tobacco", "marijuana", "sex")
Y_var  <- "birth_weight"

p <- length(X_vars)
n <- nrow(dat)

# 2. Fit the Global Outcome Model E[Y | X, C]
# We use SuperLearner to flexibly model the outcome surface
SL.lib <- c("SL.glmnet", "SL.glm", "SL.earth", "SL.xgboost", "SL.ranger")

message("Fitting global outcome regression...")
set.seed(123)
out_model_X <- outcome_model(
  Y = dat[[Y_var]], 
  X = dat[, ..X_vars], 
  C = dat[, ..C_vars], 
  mu_fitter = SL_outcome_fitter, 
  SL.lib = SL.lib,
  cvControl = list(V = 10) # 5-fold CV for SuperLearner
)

# 3. Estimate Pseudo-Outcomes
# Evaluate the outcome model to extract the pseudo-outcomes (mu_X) using Regression Adjustment
message("Estimating pseudo-outcomes...")
mu_X_obj <- estimate_ERS(
  Y = dat[[Y_var]],
  X = dat[, ..X_vars],
  C = dat[, ..C_vars],
  x_eval = dat[, ..X_vars],
  out_model = out_model_X,
  estimator = "RA",
  return_vector = TRUE # Returns just the vector of predictions
)

dat$mu_X <- mu_X_obj

# 4. Dimension Reduction via MAVE
message("Performing dimension reduction...")
form <- reformulate(X_vars, response = "mu_X")

# Fit MAVE and use cross-validation to select the optimal dimension
mave_fit <- MAVE::mave(form, data = dat, method = "meanMAVE", max.dim = p)
dhat <- MAVE::mave.dim(mave_fit)


# The cross-validation is run on dimensions of 0 1 2 3 4 
# Dimension	0 	1 	2 	3 	4 	
# CV-value	3296.32 	5.25 	5.58 	6.3 	7.45 

# Extract the orthonormal projection matrix (beta) for the selected dimension
beta <- mave_fit$dir[[dhat$dim.min]] |>
  qr() |>
  qr.Q()

d <- ncol(beta)
rownames(beta) <- X_vars

# 5. Transform Exposures into Low-Dimensional Index (Z)
message(sprintf("Selected optimal dimension: d = %d", d))
Z <- as.matrix(dat[, ..X_vars]) %*% beta
Z_vars <- paste0("Z", 1:d)
dat[, (Z_vars) := as.data.table(Z)]

# 6. Final Causal Evaluation on Z using Cross-Fitting and DR
message("Estimating final causal dose-response curve with cross-fitting (L=5)...")

if (d_assoc == 1) {
  # Standard 1D grid
  z_assoc_grid_vals <- seq(quantile(dat$Z_assoc1, 0.05), quantile(dat$Z_assoc1, 0.95), length.out = 50)
  z_assoc_eval_df <- data.frame(Z_assoc1 = z_assoc_grid_vals)
  
} else if (d_assoc == 2) {
  # 2D Mesh Grid for a Surface/Contour Plot
  z1_vals <- seq(quantile(dat$Z_assoc1, 0.05), quantile(dat$Z_assoc1, 0.95), length.out = 30)
  z2_vals <- seq(quantile(dat$Z_assoc2, 0.05), quantile(dat$Z_assoc2, 0.95), length.out = 30)
  
  # expand.grid creates every possible combination of Z1 and Z2 (900 rows)
  z_assoc_eval_df <- expand.grid(Z_assoc1 = z1_vals, Z_assoc2 = z2_vals)
  
} else {
  # For d > 2, a full mesh grid becomes computationally explosive (Curse of Dimensionality).
  # In this scenario, it is often best to simply evaluate the surface at the OBSERVED data points
  # rather than an artificial grid.
  z_assoc_eval_df <- dat[, ..Z_assoc_vars]
}

# Cross-fit the curve (or surface) along the Association dimensions
assoc_ERS <- crossfit_ERS(
  Y = dat[[Y_var]],
  X = dat[, ..Z_assoc_vars],
  C = dat[, ..C_vars],
  x_eval = z_assoc_eval_df,
  estimator = "DR",
  L = 5,
  outcome_fitter = SL_outcome_fitter,
  gps_fitter = mvn_fitter,
  optimize_bw = TRUE,
  seed = 42
)

# Estimate the final curve over the grid using 5-fold cross-fitting
# This automatically trains NN_outcome_fitter and mvn_fitter inside each fold
final_ERS <- crossfit_ERS(
  Y = dat[[Y_var]],
  X = dat[, ..Z_vars],
  C = dat[, ..C_vars],
  x_eval = z_eval_df,
  estimator = "DR",
  L = 5,
  outcome_fitter = SL_outcome_fitter,  # Neural net for smooth outcome estimation
  gps_fitter = mvn_fitter,             # Multivariate normal GPS
  optimize_bw = TRUE,                  # Dynamically calculates AMSE-optimal bandwidth
  seed = 42                            # For reproducible CV folds
)

ers_plt <- ggplot(final_ERS$results, aes(x = -Z1, y = estimate)) + 
  # 1. Subtle density marks (rug plot) using the observed data
  geom_rug(data = dat, aes(x = Z1), inherit.aes = FALSE, 
           alpha = 0.2, sides = "b", length = unit(0.05, "npc")) +
  
  # 2. Main curve and CI ribbon
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.3, fill = "steelblue") + 
  geom_line(linewidth = 1.2, color = "darkblue") + 
  
  # 3. Slide-ready theme and text sizes
  theme_bw(base_size = 18) + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text = element_text(color = "black")
  ) +
  
  # 4. Clean mathematical x-axis label
  xlab(expression(Z == 0.77*PFOS + 0.14*PFOA + 0.26*PFNA - 0.56*PFHxS )) + 
  ylab("Birthweight (g)") + 
  xlim(c(-1.8, 2.1)) +  
  ggtitle(NULL)

ers_plt
ggsave(filename = "results/jasa-initial-submission/ATL-AA/causal_ers_no_title.pdf",
       plot = ers_plt, width = 8, height = 6)


# MAVE --------------------------------------------------------------------

message("Performing Association Dimension Reduction (Standard MAVE)...")

# 1. Direct formula: Regress Y directly on X (skip pseudo-outcomes)
form_assoc <- reformulate(X_vars, response = Y_var)

# 2. Fit standard MAVE and select dimension
mave_fit_assoc <- MAVE::mave(form_assoc, data = dat, method = "meanMAVE", max.dim = p)
dhat_assoc <- MAVE::mave.dim(mave_fit_assoc)

# 3. Extract the naive association projection matrix
#beta_assoc <- mave_fit_assoc$dir[[dhat_assoc$dim.min]] |>
beta_assoc <- mave_fit_assoc$dir[[dhat$dim.min]] |>
  qr() |>
  qr.Q()

d_assoc <- ncol(beta_assoc)
rownames(beta_assoc) <- X_vars

# 4. Transform Exposures into Association Index (Z_assoc)
message(sprintf("Selected optimal association dimension: d = %d", d_assoc))
Z_assoc_mat <- as.matrix(dat[, ..X_vars]) %*% beta_assoc
Z_assoc_vars <- paste0("Z_assoc", 1:d_assoc)
dat[, (Z_assoc_vars) := as.data.table(Z_assoc_mat)]

# Association
z_assoc_grid_vals <- seq(quantile(dat$Z_assoc1, 0.01), quantile(dat$Z_assoc1, 0.95), length.out = 100)
z_assoc_eval_df <- data.frame(Z_assoc1 = z_assoc_grid_vals)

# Handle multi-dimensional grids if MAVE selected d > 1
if (d_assoc > 1) {
  for (dim_idx in 2:d_assoc) {
    z_assoc_eval_df[[paste0("Z_assoc", dim_idx)]] <- median(dat[[paste0("Z_assoc", dim_idx)]])
  }
}

# Cross-fit the curve along the Association dimension
assoc_ERS <- crossfit_ERS(
  Y = dat[[Y_var]],
  X = dat[, ..Z_assoc_vars],
  C = dat[, ..C_vars],
  x_eval = z_assoc_eval_df,
  estimator = "DR",
  L = 5,
  outcome_fitter = SL_outcome_fitter,
  gps_fitter = mvn_fitter,
  optimize_bw = TRUE,
  seed = 42
)

assoc_plt <- ggplot(assoc_ERS$results, aes(x = Z_assoc1, y = estimate)) + 
  # 1. Subtle density marks (rug plot) using the observed data
  geom_rug(data = dat, aes(x = Z_assoc1), inherit.aes = FALSE, 
           alpha = 0.2, sides = "b", length = unit(0.05, "npc")) +
  
  # 2. Main curve and CI ribbon
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.3, fill = "steelblue") + 
  geom_line(linewidth = 1.2, color = "darkblue") + 
  
  # 3. Slide-ready theme and text sizes
  theme_bw(base_size = 18) + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text = element_text(color = "black")
  ) +
  
  # 4. Clean mathematical x-axis label
  xlab(expression(Z == 0.27*PFHxS - 0.41*PFOS - 0.83*PFOA - 0.27*PFNA)) + 
  ylab("Birthweight (g)") + 
  xlim(c(-1.95, 2.1)) +  
  #ylim(c(3200, 3651)) + 
  ggtitle("Association-based Exposure Response Surface")
assoc_plt


# Beta
# Pi(Beta)
# diagonal of Pi(beta)

# Causal ERS
ggsave(filename = "results/jasa-initial-submission/ATL-AA/assocation_ers.pdf",
       plot = assoc_plt, width = 8, height = 6)

diag(Pi(beta))
diag(Pi(beta_assoc))

source("R/pCCA.R")


# Heatmap -----------------------------------------------------------------

library(ggplot2)
library(dplyr)

# 1. Create the data frame with your exact values
plot_data <- data.frame(
  Exposure = rep(c("PFOS", "PFOA", "PFNA", "PFHxS"), 2),
  Method = rep(c("CSDR (d=1)", "SDR (d=1)"), each = 4),
  Value = c(
   diag(Pi(beta)), # Causal SDR values
    diag(Pi(beta_assoc))  # Standard SDR values
  )
)

# 2. Reorder factor levels so "Causal SDR" appears on the top row, 
# and the exposures maintain their original order left-to-right.
plot_data$Method <- factor(plot_data$Method, levels = c("SDR (d=1)", "CSDR (d=1)"))
plot_data$Exposure <- factor(plot_data$Exposure, levels = c("PFOS", "PFOA", "PFNA", "PFHxS"))

# 3. Create the heatmap
heatmap <- ggplot(plot_data, aes(x = Exposure, y = Method, fill = Value)) +
  # Draw the squares with a small white border
  geom_tile(color = "white", linewidth = 1) +
  
  # Add the numeric values inside the squares (rounded to 3 decimal places)
  # Dynamically switch text color to white if the background is too dark
  geom_text(aes(label = sprintf("%.3f", Value), 
                color = Value > 0.6), 
            size = 6, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = c("black", "white")) +
  
  # Force the color scale to be exactly between 0 and 1
  scale_fill_gradient(low = "white", high = "darkblue", limits = c(0, 1)) +
  
  # Clean, slide-ready theme
  theme_minimal(base_size = 18) +
  labs(
    #title = "Subspace importance score",
    title = NULL,
    x = NULL,
    y = NULL,
    fill = "Importance\n(0 to 1)"
  ) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold", color = "black"),
    axis.text.y = element_text(face = "bold", color = "black")
  )

ggsave(filename = "results/jasa-initial-submission/ATL-AA/diag_plot_notitle.pdf",
       plot = heatmap, width = 8, height = 4)


library(ggplot2)
library(dplyr)
library(tidyr)

# 1. Define the methods and dimensions
methods_list <- c("PCA", "pCCA", "MAVE", "Oracle-MAVE", "RA-MAVE", "DR-MAVE", "PO-MAVE", "RP-MAVE")
dims_list <- paste0("Dim ", 1:10)

# Hardcode the exact mean values from the screenshot row by row
means_matrix <- matrix(c(
  0.002, 0.003, 0.013, 0.018, 0.969, 0.979, 0.007, 0.004, 0.003, 0.002, # PCA
  0.546, 0.716, 0.134, 0.134, 0.002, 0.002, 0.132, 0.122, 0.132, 0.081, # pCCA
  0.530, 0.663, 0.370, 0.022, 0.002, 0.001, 0.117, 0.110, 0.113, 0.072, # MAVE
  1.000, 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, # Oracle-MAVE
  0.880, 0.947, 0.050, 0.006, 0.000, 0.000, 0.033, 0.030, 0.031, 0.022, # RA-MAVE
  0.860, 0.928, 0.042, 0.005, 0.000, 0.001, 0.048, 0.045, 0.042, 0.029, # DR-MAVE
  0.770, 0.874, 0.044, 0.007, 0.001, 0.001, 0.088, 0.082, 0.079, 0.054, # PO-MAVE
  0.500, 0.580, 0.155, 0.163, 0.002, 0.002, 0.164, 0.171, 0.157, 0.105  # RP-MAVE
), byrow = TRUE, nrow = 8, dimnames = list(methods_list, dims_list))

# Convert the matrix into a long-format data frame for ggplot
plot_data <- as.data.frame(as.table(means_matrix))
colnames(plot_data) <- c("Method", "Dimension", "Value")

# 2. Reorder factor levels 
# rev() ensures "PCA" appears on the top row instead of the bottom
plot_data$Method <- factor(plot_data$Method, levels = rev(methods_list))
plot_data$Dimension <- factor(plot_data$Dimension, levels = dims_list)

# 3. Create the heatmap
heatmap <- ggplot(plot_data, aes(x = Dimension, y = Method, fill = Value)) +
  # Draw the squares with a small white border
  geom_tile(color = "white", linewidth = 1) +
  
  # Add the numeric values inside the squares (rounded to 3 decimal places)
  # Dynamically switch text color to white if the background is too dark
  geom_text(aes(label = sprintf("%.3f", Value), 
                color = Value > 0.6), 
            size = 3.5, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = c("black", "white")) +
  
  # Force the color scale to be exactly between 0 and 1
  scale_fill_gradient(low = "white", high = "darkblue", limits = c(0, 1)) +
  
  # Clean, slide-ready theme
  theme_minimal(base_size = 18) +
  labs(
    title = NULL,
    x = NULL,
    y = NULL,
    fill = "Importance\n(0 to 1)"
  ) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold"),
    # Angled x-axis text to prevent overlap with 10 columns
    axis.text.x = element_text(face = "bold", color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(face = "bold", color = "black")
  )

# Display the plot
print(heatmap)

# 4. Save the plot
ggsave(filename = "results/jasa-initial-submission/ATL-AA/diag_plot_n=500.pdf",
        plot = heatmap, width = 12, height = 5)
library(ggplot2)
library(dplyr)
library(tidyr)

# 1. Define the methods and dimensions
methods_list <- c("PCA", "pCCA", "MAVE", "Oracle-MAVE", "RA-MAVE", "DR-MAVE", "PO-MAVE", "RP-MAVE")
dims_list <- paste0("Dim ", 1:10)

# Hardcode the exact mean values from the n=5000 screenshot row by row
means_matrix <- matrix(c(
  0.002, 0.003, 0.011, 0.015, 0.972, 0.982, 0.007, 0.004, 0.003, 0.002, # PCA
  0.819, 0.884, 0.053, 0.053, 0.001, 0.001, 0.053, 0.053, 0.053, 0.031, # pCCA
  0.363, 0.782, 0.764, 0.010, 0.002, 0.000, 0.023, 0.021, 0.022, 0.013, # MAVE
  1.000, 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, # Oracle-MAVE
  0.999, 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, # RA-MAVE
  0.996, 0.998, 0.002, 0.000, 0.000, 0.000, 0.001, 0.001, 0.001, 0.001, # DR-MAVE
  0.945, 0.978, 0.001, 0.001, 0.000, 0.000, 0.018, 0.021, 0.022, 0.013, # PO-MAVE
  0.827, 0.878, 0.050, 0.053, 0.001, 0.001, 0.052, 0.053, 0.052, 0.033  # RP-MAVE
), byrow = TRUE, nrow = 8, dimnames = list(methods_list, dims_list))

# Convert the matrix into a long-format data frame for ggplot
plot_data <- as.data.frame(as.table(means_matrix))
colnames(plot_data) <- c("Method", "Dimension", "Value")

# 2. Reorder factor levels 
# rev() ensures "PCA" appears on the top row instead of the bottom
plot_data$Method <- factor(plot_data$Method, levels = rev(methods_list))
plot_data$Dimension <- factor(plot_data$Dimension, levels = dims_list)

# 3. Create the heatmap
heatmap <- ggplot(plot_data, aes(x = Dimension, y = Method, fill = Value)) +
  # Draw the squares with a small white border
  geom_tile(color = "white", linewidth = 1) +
  
  # Add the numeric values inside the squares (rounded to 3 decimal places)
  # Dynamically switch text color to white if the background is too dark
  geom_text(aes(label = sprintf("%.3f", Value), 
                color = Value > 0.6), 
            size = 3.5, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = c("black", "white")) +
  
  # Force the color scale to be exactly between 0 and 1
  scale_fill_gradient(low = "white", high = "darkblue", limits = c(0, 1)) +
  
  # Clean, slide-ready theme
  theme_minimal(base_size = 18) +
  labs(
    title = NULL,
    x = NULL,
    y = NULL,
    fill = "Importance\n(0 to 1)"
  ) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold"),
    # Angled x-axis text to prevent overlap with 10 columns
    axis.text.x = element_text(face = "bold", color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(face = "bold", color = "black")
  )

# Display the plot
print(heatmap)

# 4. Save the plot
ggsave(filename = "results/jasa-initial-submission/ATL-AA/diag_plot_n=5000.pdf",
        plot = heatmap, width = 12, height = 5)
