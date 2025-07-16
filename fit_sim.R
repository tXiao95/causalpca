library(here)

source(here("R", "gcomp.R"))
source(here("R", "kernel_smooth.R"))
source(here("R", "csPCA.R"))
source(here("R", "doPCA.R"))
source(here("R", "misc.R"))
source(here("R", "sim_data.R"))
source(here("R", "population_obj.R"))

# Main code ---------------------------------------------------------------
main <- function(){
  set.seed(TASK_ID)
  # sim        <- sim_linear(n = 200, 
  #                          Sigma_e = diag(c(9, 1)), 
  #                          beta = c(1,0.1)
  #                          )
  
  sim   <- sim_unif_e(n=200)
  
  # # Empirical covariance matrices. Cov(X,Y) is p x 1 and Sigma_XX is p x p.
  Sigma_XY <- cov(sim$X, sim$Y) |> as.numeric()
  Sigma_XX <- cov(sim$X)
  
  # # Associational PCA
  omega.pca <- prcomp(sim$X, center = TRUE, scale. = FALSE)$rotation[,1]
  omega.pls <- unitvec(Sigma_XY)                             # PLS
  omega.cca <- unitvec(solve(Sigma_XX, Sigma_XY))            # CCA
  
  # Associational PCA w/ confounders
  # fit.lm     <- lm(sim$Y ~ sim$X + sim$C + 0)
  # omega.ppls <- 
  # omega.pcca <- solve(Sigma_XX) %*% omega.ppls |> unitvec()
  
  # Causal PCA
  cs   <- csPCA(sim$Y, sim$X, sim$C, maxit = 1000)
  do   <- doPCA(sim$Y, sim$X, sim$C, maxit = 1000)
  do.s <- doPCA(sim$Y, sim$X, sim$C, maxit = 1000, scaled = TRUE)

  # Organize estimates
  est <- rbind(cs$omega, do$omega, do.s$omega,
               omega.pca, omega.pls, omega.cca) |> data.frame()
  
  colnames(est) <- paste0("X", 1:ncol(sim$X))
  est[est$X1 < 0, ] <- -est[est$X1 < 0, ]
  
  est$method    <- c("cs", "do", "do_scaled", "PCA", "PLS", "CCA")
  est$ID        <- TASK_ID
  
  # Final object
  obj <- list(est = est, cs = cs, do = do, do.s = do.s, data = sim)
  
  saveRDS(obj, paste0("outputs_unif/", TASK_ID, ".rds"))
}

if(interactive()){
  TASK_ID    <- 1
} else{
  TASK_ID    <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
}

main()
