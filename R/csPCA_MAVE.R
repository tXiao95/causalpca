source("R/gcomp.R")
library(MAVE)

csPCA <- function(Y, X, C, ...){
  n    <- nrow(X); p <- ncol(X); q <- ncol(C)
  # Replace mu_X with any estimate of causal mean mu(X). Here we use gcomp or regression adjustment. 
  mu_X <- gcomp(Y, X, C)
  mave <- MAVE::mave(mu_X ~ as.matrix(X), method = "meanMAVE", ...)
  
  return(mave)
}