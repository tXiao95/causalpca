library(data.table)
library(here)

loadings <- readRDS(here("outputs/cytokines/beta_varimax_loadings.rds"))

L <- loadings$`MC__TNF‑α_c`

varimax(L$beta_final)
