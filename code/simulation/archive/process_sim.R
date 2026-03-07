library(data.table)
library(ggplot2)
library(here)

# Process -----------------------------------------------------------------
#dir      <- paste0("enar-student-paper2_", corr, "-rho", "_p", nvar)
dir      <- paste0("enar-student-paper2_highZ1")
sim_path <- here("outputs", "simulation", dir)
dt       <- lapply(list.files(sim_path, full.names = TRUE), fread) |>
                   rbindlist(fill = TRUE)

dt[, method := ifelse(method == "CS", "csPCA (L=1)",
                      ifelse(method == "CS-CF", "csPCA (L=5)", method))]

dt.summary <- dt[, .(mean = mean(value), sd = sd(value),
       time_avg = mean(time), time_sd = sd(time)), .(n, error_type, method)]

dt[, method := ifelse(method == "Truth", "csPCA (Oracle)", method)]
dt[, n := paste0("n = ", n)]
dt[, n := factor(n, levels = c("n = 100", 
                               "n = 200", 
                               "n = 400", 
                               "n = 800", 
                               "n = 1600"))]

# Boxplot of errors -------------------------------------------------------
ggplot(dt[n %in% c("n = 100","n = 400","n = 1600") & error_type == "Frobenius"], aes(method, value)) + 
  geom_boxplot(aes(col = method)) + 
  #facet_grid(error_type~n, scales = "free") + 
  facet_wrap(~n) + 
  theme_bw(base_size = 16) + 
  ylab("Frobenius error") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) 
ggsave("results/simulation-boxplot.png", width = 11, height = 8.5)
ggsave("results/simulation-boxplot.pdf", width = 11, height = 5)

# Table of Errors (Frobenius) ---------------------------------------------------------
dt.summary[, value := paste0(round(mean, 4), " (", round(sd, 4), ")")]
dt.wide <- dt.summary[error_type == "2-norm"] |> 
  dcast(n ~ method, value.var = "value")
dt.wide <- dt.summary[error_type == "Frobenius"] |> 
  dcast(n ~ method, value.var = "value")

dt.wide[, .(n, PCA, pCCA, SDR, `Truth`, `csPCA (L=1)`, `csPCA (L=5)`)]

# Table of Runtimes (minutes) -------------------------------------------------------
dt.summary[, value := paste0(round(time_avg / 60, 4), " (", round(time_sd / 60, 4), ")")]
dt.wide <- dt.summary[error_type == "Frobenius"] |> 
  dcast(n ~ method, value.var = "value")

dt.wide[, .(n, PCA, pCCA, SDR, `csPCA (L=1)`, `csPCA (L=5)`)]

# Sqrt(n) Error ()-----------------------------------------------------------
ggplot(dt.summary[method != "Truth"], aes(n^{-1/2}, mean)) + 
  geom_point() + 
  geom_line(aes(group = method, col = method)) +
  facet_wrap(~error_type, scales = "free") + 
  xlab("Sample size (n)") + 
  ylab("Error") + 
  theme_bw() + ggtitle("Error") + 
  #ylim(c(0, 2))  + 
  xlim(c(0, .1)) + 
  #geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0)
ggsave("results/simulation-avg-error.png")

# MSE of mu(X) ------------------------------------------------------------
dt.X <- dt[method %in% c("SDR", "csPCA (L=1)", "csPCA (L=5)"),
           .(value=paste0(mean(mu_X_mse)," (", sd(mu_X_mse), ")")), 
                                                                 .(n, method)] |>
  dcast(n ~ method, value.var = "value")

dt.X

# MSE of mu(Z) ------------------------------------------------------------
dt.Z <- dt[method %in% c("Truth", "SDR", "csPCA (L=1)", "csPCA (L=5)"),
           .(value=paste0(mean(mu_Z_mse)," (", sd(mu_Z_mse), ")")), 
                                                                 .(n, method)] |>
  dcast(n ~ method, value.var = "value")

dt.Z

# Dimension selection -----------------------------------------------------
nd <- length(unique(dt$ID))*2
dt[method %in% c("SDR", "csPCA (L=1)", "csPCA (L=5)", "Truth") & (n ==100), dhat, .(n, method)] |>
  dcast(method ~ dhat, fun.aggregate = function(x) length(x) / nd)
dt[method %in% c("SDR", "csPCA (L=1)", "csPCA (L=5)", "Truth") & (n ==200), dhat, .(n, method)] |>
  dcast(method ~ dhat, fun.aggregate = function(x) length(x) / nd)
dt[method %in% c("SDR", "csPCA (L=1)", "csPCA (L=5)", "Truth") & (n ==400), dhat, .(n, method)] |>
  dcast(method ~ dhat, fun.aggregate = function(x) length(x) / nd)
dt[method %in% c("SDR", "csPCA (L=1)", "csPCA (L=5)", "Truth") & (n ==800), dhat, .(n, method)] |>
  dcast(method ~ dhat, fun.aggregate = function(x) length(x) /nd)
dt[method %in% c("SDR", "csPCA (L=1)", "csPCA (L=5)", "Truth") & (n ==1600), dhat, .(n, method)] |>
  dcast(method ~ dhat, fun.aggregate = function(x) length(x) / nd)

# MSE of mu(Z) ------------------------------------------------------------

dt[, .(mean(mu_Z_mse, na.rm = TRUE), sd(mu_Z_mse, na.rm = TRUE)), .(n, method)]
dt.xmse <- dt[method %in% c("CS", "CS-CF", "SDR", "Truth"), 
   .(mse=mean(mu_X_mse, na.rm = TRUE), sd=sd(mu_X_mse, na.rm = TRUE)), .(n, method)]

dt.zmse <- dt[method %in% c("csPCA (L=1)", "csPCA (L=5)", "SDR", "Truth"), 
   .(mse=mean(mu_Z_mse, na.rm = TRUE), sd=sd(mu_Z_mse, na.rm = TRUE)), .(n, method)]

dt.zmse[, value := paste0(round(mse, 4), " (", round(sd, 4), ")")]

a[, lower := mse - 1.96*sd]
a[, upper := mse + 1.96*sd]

ggplot(a, aes(n, mse)) + 
  geom_point()+ 
  #geom_errorbar(aes(ymin = lower, ymax = upper)) + 
  geom_line(aes(group = method, col = method)) + 
  #ylim(c(0,0.3)) + 
  geom_hline(yintercept = 0) + 
  theme_bw() + 
  xlim(c(0,1600)) + 
  geom_vline(xintercept = 0) + 
  ylab("MSE of mu(Z)") + 
  xlab("Sample size (n)")  
  #facet_wrap(~method)

