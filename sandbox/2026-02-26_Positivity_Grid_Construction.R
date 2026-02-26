# Date: 2026/02/26

# Observe differences in how to make a positivity grid for observing effects
# The main point is we want to avoid positivity violations. These lead to 
# potentially disastrous extrapolation in the real data analysis.

library(here)
source(here("R/get_positivity_grid.R"))

set.seed(42)
mean_vec <- c(0, 0)
cov_mat <- matrix(c(1, 0.9, 0.9, 1), nrow = 2) # High correlation "football" shape
data <- MASS::mvrnorm(n = 500, mu = mean_vec, Sigma = cov_mat)
colnames(data) <- c("Var1", "Var2")

# Generate the tight grid
kde_grid    <- get_kde_positivity_grid(data, percentile = 0, n_points = 30, threshold = 0.1)
convex_grid <- get_convex_positivity_grid(data, percentile = 0, n_points = 30)


pdf(here("outputs/experiments/Positivity_Grid_Construction.pdf"), width = 8, height = 6)
# --- Visualization ---
# 1. First, create the empty plot with the original data
plot(data, col = rgb(0, 0, 0, 0.5), pch = 16, 
     main = "KDE Density (Red) vs. Convex Hull (Blue)",
     xlab = "Variable 1", ylab = "Variable 2")

# 2. Add the Convex Hull grid in Blue
# (Note: Fixed the variable name typo from 'convex_grid_grid')
points(convex_grid, col = "blue", pch = 15, cex = 0.8)

# 3. Add the Tight KDE grid in Red
# Since this is a subset of the convex grid, these points will 
# overlap the blue ones, effectively 'changing' their color to red.
points(kde_grid, col = "red", pch = 14, cex = 0.8)

# 4. Add a legend to keep track of which is which
legend("topleft", 
       legend = c("Original Data", "Convex Hull Grid", "Tight KDE Grid"), 
       col = c("gray", "blue", "red"), 
       pch = c(16, 15, 14),
       bty = "n") # Removes the box around the legend
dev.off()