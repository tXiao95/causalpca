library(here)
library(data.table)
library(ggplot2)
library(patchwork)
source(here("R/misc.R"))

# read in objects labeled "bw_obj_N=XXX.rds" where XXX is the sample size as an integer
filepaths <- list.files(here("results", "birthweight"), "bw_obj_*", full.names = TRUE)

N   <- as.numeric(sub(".*N=([0-9]+)\\.rds$", "\\1", filepaths))
idx <- order(N)

N <- N[idx]
filepaths <- filepaths[idx]

results <- lapply(filepaths, function(path){
  print(path)
  N   <- as.numeric(sub(".*N=([0-9]+)\\.rds$", "\\1", path))
  obj <- readRDS(path)
  
  out_model_X_time <- obj$out_model_X_time[3]
  gcomp_time <- obj$gcomp_time[3]
  mave_time <- obj$mave_time[3]
  dhat_time <- obj$dhat_time[3]
  out_model_Z_time <- obj$out_model_Z_time[3]
  
  
  P <- obj$P
  d <- obj$dhat$dim.min
  runtimes <- data.table(N = N, 
                         type = c("Outcome X",
                                  "gcomp",
                                  "MAVE",
                                  "dhat",
                                  "Outcome Z"), 
                         time = c(out_model_X_time,
                                  gcomp_time,
                                  mave_time,
                                  dhat_time,
                                  out_model_Z_time))
  
  runtimes[, time_min := as.numeric(time) / 60]
  runtimes[, time_hour := time_min / 60]
  
  return(list(P = P, d = d, runtimes = runtimes))
})

names(results) <- N

tab <- rbindlist(lapply(results, function(x) x$runtimes))
P_list <-lapply(results, function(x) x$P)
d_list <- sapply(results, function(x) x$d)

dists <- sapply(seq_len(length(P_list) - 1),
                function(i) Matrix::norm(P_list[[i]] - P_list[[i + 1]], type = "2"))
d_tab <- data.table(d = d_list, N = N, dist = c(NA, dists))

runtime_plot <- ggplot(tab, aes(N, time_min)) +
  geom_point() + 
  geom_line() + 
  facet_wrap(~type, scales = "free") + 
  #scale_x_log10() + 
  #scale_y_log10() + 
  ylab("Minutes") + 
  xlab("N (linear scale)") + 
  #geom_smooth(method = "loess", alpha = 0.2) + 
  ggtitle("Runtimes (min) for different components of causal SDR")

d_plot <- ggplot(d_tab, aes(N, d)) + 
  geom_point() + 
  geom_line() + 
  scale_x_log10() +
  scale_y_continuous(breaks = seq(1, max(d_tab$d), by = 1)) + 
  #theme_minimal() + 
  ylab("dhat") + 
  ggtitle("Structural dimension selection")

dist_plot <- ggplot(d_tab[2:nrow(d_tab),], aes(N, dists)) + 
  geom_point() + 
  geom_line() + 
  scale_x_log10() +
  ylab("2-norm") + 
  xlab("N (log scale)") + 
  coord_cartesian(ylim=c(0,1)) + 
  ggtitle("Distances between Projection Matrices")

(d_plot | dist_plot) / runtime_plot + 
  plot_layout(heights = c(1,2))

ggsave(filename = here("results", "birthweight", "BW_runtimes.pdf"), width = 11, height = 8.5)

