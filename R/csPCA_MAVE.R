source("R/gcomp.R")
library(MAVE)

csPCA <- function(Y, X, C, ..., gcomp_args = list(), mave_args = list()) {
  # gcomp call
  mu_X <- do.call(gcomp, c(list(Y = Y, X = X, C = C), gcomp_args))
  
  # mave call
  mave <- do.call(MAVE::mave,
                  c(list(formula = mu_X ~ as.matrix(X), method = "meanMAVE"),
                    mave_args))
  
  return(mave)
}
