plot_ers_pointwise_coverage <- function(dt_ers, dt_true_z) {
  if (is.null(dt_ers)) return(NULL)
  
  dt_main <- copy(dt_ers[method %in% main_methods])
  dt_main[, method := droplevels(method)]
  
  # Pointwise coverage at each of the 100 evaluation points
  dt_pointwise <- dt_main[, .(
    coverage = mean(mu_true >= ci_lower & mu_true <= ci_upper, na.rm = TRUE)
  ), by = .(eval_id, method, n_str)]
  
  # Merge true subspace coordinates
  dt_pointwise <- merge(dt_pointwise, dt_true_z, by = "eval_id", all.x = TRUE)
  
  # Keep selected methods / sample sizes
  dt_pointwise <- dt_pointwise[
    n_str %in% c("n = 500", "n = 5000") &
      method %in% c("Oracle-MAVE", "RA-MAVE", "Full_X", "MAVE", "pCCA")
  ]
  
  # --- PLOTMATH FIXES ---
  dt_pointwise[, n_str := fifelse(n_str == "n = 500", "n == 500", "n == 5000")]
  dt_pointwise[, n_str := factor(n_str, levels = c("n == 500", "n == 5000"))]
  
  dt_pointwise[, method_label := fifelse(
    method == "Oracle-MAVE", "plain('Oracle')",
    fifelse(
      method == "RA-MAVE", "plain('Causal SDR')",
      fifelse(
        method == "Full_X", "plain('Full')~X", 
        paste0("plain('", as.character(method), "')") 
      )
    )
  )]
  
  dt_pointwise[, method_label := factor(
    method_label,
    levels = c("plain('Oracle')", "plain('Causal SDR')", "plain('Full')~X", "plain('MAVE')", "plain('pCCA')")
  )]
  
  # Average coverage across the 100 evaluation points within each panel
  panel_summary <- dt_pointwise[, .(
    mean_coverage = median(coverage, na.rm = TRUE)
  ), by = .(n_str, method_label)]
  
  # Create the label string
  panel_summary[, label := sprintf("Median = %.3f", mean_coverage)]
  
  # Find the horizontal midpoint for the labels
  x_mid <- mean(range(dt_pointwise$Z1_true, na.rm = TRUE))
  
  p <- ggplot(dt_pointwise, aes(x = Z1_true, y = Z2_true)) +
    geom_point(aes(color = coverage), size = 2, alpha = 0.8) +
    geom_label(
      data = panel_summary,
      aes(x = x_mid, y = -Inf, label = label), 
      inherit.aes = FALSE,
      hjust = 0.5, 
      vjust = 0,   
      size = 3, 
      fontface = "bold",
      label.size = 0.2,
      label.padding = unit(0.2, "lines"),
      label.r = unit(0, "lines"), 
      fill = "white"
    ) +
    facet_grid(
      n_str ~ method_label, 
      labeller = label_parsed
    ) +
    scale_color_gradientn(
      colours = c("firebrick", "gray90", "dodgerblue"),
      values = c(0, 0.95, 1),
      limits = c(0, 1),
      breaks = c(0, 0.50, 0.80, 0.95),
      labels = c("0.00", "0.50", "0.80", "0.95"),
      name = "Coverage",
      guide = guide_colorbar(
        title.position = "left", 
        title.vjust = 0.8,       
        title.hjust = 1
      )
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.25, 0.05))) + 
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
    labs(
      x = expression(True~Z[1]),
      y = expression(True~Z[2])
    ) +
    theme_bw(base_size = 9) + 
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "gray95"),
      # Explicitly lock all major titles and labels to size 11
      strip.text = element_text(face = "bold", size = 11), 
      axis.title = element_text(face = "bold", size = 11), 
      legend.title = element_text(face = "bold", size = 11),
      # Optional: explicitly set axis ticks and legend text slightly smaller so the titles stand out
      axis.text = element_text(size = 9),
      legend.text = element_text(size = 9),
      legend.position = "bottom",
      legend.key.height = unit(0.4, "cm"), 
      legend.key.width = unit(1.2, "cm"),
      legend.margin = margin(t = -5, b = 0) 
    )
  
  return(p)
}

plot <- plot_ers_pointwise_coverage(dt_ers, dt_true_z)

ggsave(
  here("results/jasa-initial-submission/main_paper_final_results_nnet/interaction/plot_pw_coverage_main_paper.pdf"), 
  width = 7.2, 
  height = 4.2, 
  dpi = 300
)
ggsave(
  here("results/jasa-initial-submission/main_paper_final_results_nonsparse_nnet_128x64x32/interaction/plot_pw_coverage_main_paper.pdf"), 
  width = 7.2, 
  height = 4.2, 
  dpi = 300
)
