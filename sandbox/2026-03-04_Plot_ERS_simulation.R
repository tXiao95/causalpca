# 1. Define the true causal mean function
h_Z_fun <- function(Z1, Z2) {
  return( 3 * tanh(Z1) + 8 * (pnorm(2 * Z2) - 0.5) )
}

# 2. Create a grid of values for the structural dimensions Z1 and Z2
# Ranging from -3 to 3 covers the area where both functions level off
z1_vals <- seq(-3, 3, length.out = 50)
z2_vals <- seq(-3, 3, length.out = 50)

# 3. Evaluate the function across the grid
h_Z_grid <- outer(z1_vals, z2_vals, Vectorize(h_Z_fun))

# 4. Generate the static 3D surface plot
persp(x = z1_vals, 
      y = z2_vals, 
      z = h_Z_grid,
      theta = 45, phi = 30, expand = 0.6,
      col = "cadetblue1", ltheta = 120, shade = 0.75,
      ticktype = "detailed",
      xlab = "Z1 (Causal Direction 1)",
      ylab = "Z2 (Causal Direction 2)",
      zlab = "Causal Effect h(Z)",
      main = "True Causal Response Surface")
