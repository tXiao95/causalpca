library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

# 1. Define the original methods
methods_raw <- c("PCA", "pCCA", "MAVE", "Oracle-MAVE", "RA-MAVE", "DR-MAVE", "PO-MAVE", "RP-MAVE")
methods_raw <- c("PCA", "pCCA", "MAVE", "Oracle-MAVE", "RA-MAVE", "DR-MAVE", "PO-MAVE")

# Clean labels: Remove "-MAVE" but keep "MAVE" as is
methods_clean <- ifelse(methods_raw == "MAVE", "MAVE", str_replace(methods_raw, "-MAVE", ""))

# Create mathematical labels for X1 to X10
dims_list <- paste0("X", 1:10)
dims_math <- paste0("X[", 1:10, "]") 

# 2. Hardcode the exact Means and SDs from the n=500 text
means_vec <- c(
  0.002, 0.003, 0.013, 0.018, 0.969, 0.979, 0.007, 0.004, 0.003, 0.002, # PCA
  0.546, 0.716, 0.134, 0.134, 0.002, 0.002, 0.132, 0.122, 0.132, 0.081, # pCCA
  0.530, 0.663, 0.370, 0.022, 0.002, 0.001, 0.117, 0.110, 0.113, 0.072, # MAVE
  1.000, 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, # Oracle
  0.880, 0.947, 0.050, 0.006, 0.000, 0.000, 0.033, 0.030, 0.031, 0.022, # RA
  0.860, 0.928, 0.042, 0.005, 0.000, 0.001, 0.048, 0.045, 0.042, 0.029, # DR
  0.770, 0.874, 0.044, 0.007, 0.001, 0.001, 0.088, 0.082, 0.079, 0.054, # PO
  0.500, 0.580, 0.155, 0.163, 0.002, 0.002, 0.164, 0.171, 0.157, 0.105  # RP
)

sds_vec <- c(
  0.001, 0.001, 0.011, 0.012, 0.013, 0.005, 0.001, 0.001, 0.001, 0.001, # PCA
  0.131, 0.089, 0.149, 0.152, 0.003, 0.003, 0.150, 0.133, 0.143, 0.100, # pCCA
  0.282, 0.247, 0.216, 0.023, 0.002, 0.002, 0.106, 0.095, 0.107, 0.070, # MAVE
  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, # Oracle
  0.209, 0.134, 0.136, 0.019, 0.001, 0.001, 0.084, 0.066, 0.070, 0.056, # RA
  0.202, 0.125, 0.117, 0.016, 0.001, 0.001, 0.086, 0.083, 0.071, 0.055, # DR
  0.277, 0.207, 0.121, 0.023, 0.002, 0.003, 0.156, 0.134, 0.129, 0.105, # PO
  0.232, 0.227, 0.145, 0.150, 0.004, 0.003, 0.147, 0.153, 0.132, 0.102  # RP
)

means_vec <- c(
  # PCA
  0.307, 0.323, 0.094, 0.063, 0.176, 0.164, 0.183, 0.224, 0.248, 0.219,
  # pCCA
  0.449, 0.634, 0.162, 0.165, 0.002, 0.002, 0.167, 0.154, 0.153, 0.113,
  # MAVE
  0.135, 0.469, 0.388, 0.243, 0.001, 0.001, 0.118, 0.267, 0.131, 0.247,
  # Oracle
  0.249, 0.250, 0.250, 0.250, 0.000, 0.000, 0.250, 0.250, 0.250, 0.250,
  # RA
  0.266, 0.263, 0.194, 0.218, 0.000, 0.000, 0.264, 0.261, 0.262, 0.271,
  # DR
  0.236, 0.258, 0.186, 0.218, 0.000, 0.000, 0.279, 0.267, 0.281, 0.273,
  # PO
  0.245, 0.270, 0.169, 0.219, 0.001, 0.001, 0.263, 0.284, 0.267, 0.282
)

sds_vec <- c(
  # PCA
  0.021, 0.018, 0.023, 0.019, 0.008, 0.004, 0.003, 0.005, 0.010, 0.013,
  # pCCA
  0.099, 0.096, 0.171, 0.176, 0.004, 0.004, 0.176, 0.167, 0.161, 0.135,
  # MAVE
  0.115, 0.168, 0.179, 0.056, 0.002, 0.002, 0.105, 0.135, 0.123, 0.107,
  # Oracle
  0.013, 0.005, 0.010, 0.003, 0.000, 0.000, 0.020, 0.007, 0.020, 0.006,
  # RA
  0.051, 0.022, 0.042, 0.015, 0.000, 0.000, 0.063, 0.028, 0.064, 0.027,
  # DR
  0.054, 0.022, 0.041, 0.016, 0.000, 0.000, 0.072, 0.033, 0.069, 0.029,
  # PO
  0.097, 0.046, 0.067, 0.030, 0.003, 0.003, 0.125, 0.080, 0.119, 0.060
)

# 3. Build the DataFrame
plot_data <- data.frame(
  Method = rep(methods_clean, each = 10),
  #Dimension = rep(dims_list, times = 8),
  Dimension = rep(dims_list, times = 7),
  Mean = means_vec,
  SD = sds_vec
)

# Combine Mean and SD into a 2-line string rounded to 2 decimal places
plot_data <- plot_data %>%
  mutate(
    Label = sprintf("%.2f\n(%.2f)", Mean, SD),
    # Lock the factor orders
    Method = factor(Method, levels = rev(methods_clean)),
    Dimension = factor(Dimension, levels = dims_list)
  )

# 4. Create the Aesthetic Heatmap
heatmap <- ggplot(plot_data, aes(x = Dimension, y = Method, fill = Mean)) +
  geom_tile(color = "white", linewidth = 1.2) +
  
  # REMOVED fontface="bold" here so the numbers are easier to read
  geom_text(aes(label = Label, color = Mean > 0.49), 
            size = 3.5, show.legend = FALSE) +
  scale_color_manual(values = c("black", "white")) +
  
  # A sleek, modern blue gradient scale
  scale_fill_gradient(low = "#F0F8FF", high = "#00204D", limits = c(0, 1),
                      breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
  
  # Convert X1...X10 into proper subscript math labels (X_1...X_10)
  scale_x_discrete(labels = parse(text = dims_math)) +
  
  theme_minimal(base_size = 10) +
  labs(
    x = NULL,
    y = NULL,
    fill = "Subspace Importance Score"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(face = "bold", color = "black", size = 12),
    axis.text.y = element_text(face = "bold", color = "black", size = 12),
    axis.ticks = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", vjust = 0.8),
    legend.key.width = unit(1.5, "cm"),
    legend.key.height = unit(0.4, "cm")
  )

print(heatmap)

# Save the updated un-bolded version
ggsave("heatmap_n500.pdf", plot = heatmap, width = 7.2, height = 5, dpi = 300)
ggsave("heatmap_n1000_nonsparse.pdf", plot = heatmap, width = 7.2, height = 5, dpi = 300)
