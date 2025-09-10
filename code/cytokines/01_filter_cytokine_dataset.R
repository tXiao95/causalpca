library(data.table)
library(forcats)
library(ggplot2)
library(here)
library(naniar)

# Filepaths ---------------------------------------------------------------
processed_data_path <- here("data/cytokines/processed/")
full_data_path <- here(processed_data_path, "full(n=106)_cytokine_mpc_log-z-transform_w_covariates.csv")

# Collapse categories --------------------------------------------------------------------

dt <- fread(full_data_path, stringsAsFactors = TRUE, na.strings = c("", "NA"))
dt[, smoke := fct_collapse(
  smoke,
  "Former or Current smoker" = c("Current smoker", "Former smoker")
)]

dt[, drugs := fct_collapse(
  drugs,
  "Former or Current user" = c("Current User", "Former User")
)]

naniar::vis_miss(dt[,40:ncol(dt)])
ggsave(here("outputs/cytokines/figures/Missingness Pattern of Cytokine Dataset.png"))


# These columns have too many missing (>=4% of the 106 observations) or are not informative (gest_multiple)
drop_cols <- c("marital", "mat_race_eth", "foreign", "gest_multiple", "smoke")

dt <- dt[, !..drop_cols]
dt <- dt[complete.cases(dt)] 

fwrite(dt, here(processed_data_path, "reduced(n=100)_cytokine_mpc_log-z-transform_w_covariates.csv"))

# Make datasets -----------------------------------------------------------

# Eventually, we want to do a whole multidimensional treatment + mediation analysis for
# the maternal -->> placental -->> cord blood pathway. For now, create three datasets for

# 1. Maternal -->> placental
# 2. Maternal -->> cord blood
# 3. Placental -->> cord blood

m_cols <- names(dt)[grep("_m",names(dt))]
p_cols <- names(dt)[grep("_p",names(dt))]
c_cols <- names(dt)[grep("_c",names(dt))]

p_outcome_cols <- c("TNF‑α_p", "IL-6_p")
c_outcome_cols <- c("TNF‑α_c", "IL-6_c")

p_left <- setdiff(p_cols, p_outcome_cols)
c_left <- setdiff(c_cols, c_outcome_cols)

# Filter out extra columns. Only keep outcome columns. 
MP <- dt[, !..c_cols][, !..p_left]
MC <- dt[, !..p_cols][, !..c_left]
PC <- dt[, !..m_cols][, !..c_left]

# Write to disk
fwrite(MP, here(processed_data_path, "analysis_maternal-to-placental.csv"))
fwrite(MC, here(processed_data_path, "analysis_maternal-to-cord.csv"))
fwrite(PC, here(processed_data_path, "analysis_placental-to-cord.csv"))