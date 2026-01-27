library(data.table)
library(ggplot2)
library(here)
library(stringr)

# Filepaths ---------------------------------------------------------------
results_path        <- here("outputs/cytokines/figures/")
cytokine_data_path  <- here("data/cytokines/raw/Cytokines/WOC cytokines_cleaned_06.24.2024.csv")
covariate_data_path <- here("data/cytokines/raw/Covariate and stress/WOC covar and stress data_05022024.csv")
final_data_path     <- here("data/cytokines/processed/")

# Functions ---------------------------------------------------------------
# Reads in cytokines in long format: easier for plotting and preprocessing
# type: "m" is maternal, "p" is placental, "c" is cord blood
reshape_cytokines <- function(type = c("m","p","c"),
                              file = cytokine_data_path) {
  type <- match.arg(type)
  
  DT <- fread(file)
  DT <- DT[, crh_m := NULL] # Not sure what this is...variable not in codebook
  
  # --- helper: strip suffix & build a matching key (case/sep-insensitive) ---
  strip_suffix <- function(x, type) {
    x <- sub(paste0("_", type, "$"), "", x)
    x <- sub(paste0("\\.", type, "\\.LOD$"), "", x)
    x
  }
  keyify <- function(x) {
    x |>
      str_to_lower() |>
      str_replace_all("[^a-z0-9]", "")  # keep only letters/digits
  }
  
  # --- canonical base-name mapping for all 20 cytokines ---
  canon <- c(
    il28a   = "IL-28A",  il28   = "IL-28A",
    il6     = "IL-6",
    il8     = "IL-8",
    il7     = "IL-7",
    cxcl10  = "CXCL10",
    il10    = "IL-10",
    ccl2    = "CCL2",
    il1beta = "IL-1Beta",
    ifngamma= "IFN-Gamma",
    ccl3    = "CCL3",
    ccl22   = "CCL22",
    ccl4    = "CCL4",
    il4     = "IL-4",
    il2     = "IL-2",
    cxcl9   = "CXCL9",
    ifnalpha= "IFN-Alpha",
    tnfalpha= "TNF-Alpha",
    tnfri   = "TNF RI",
    il12    = "IL-12",
    ccl17   = "CCL17"
  )
  
  # --- select columns for this type ---
  val_pat <- paste0("_", type, "$")
  lod_pat <- paste0("\\.", type, "\\.LOD$")
  
  keep_val <- grep(val_pat, names(DT), value = TRUE)
  keep_lod <- grep(lod_pat, names(DT), value = TRUE)
  
  # always keep ppt_id
  sel <- c("ppt_id", keep_val, keep_lod)
  sel <- sel[sel %in% names(DT)]
  DT  <- DT[, ..sel]
  
  # --- normalize names so bases match across value & lod ---
  # build canonical base names for value columns
  bases_val <- strip_suffix(keep_val, type)
  keys_val  <- keyify(bases_val)
  canon_val <- unname(canon[keys_val])
  
  # build canonical base names for lod columns (if any)
  if (length(keep_lod)) {
    bases_lod <- strip_suffix(keep_lod, type)
    keys_lod  <- keyify(bases_lod)
    canon_lod <- unname(canon[keys_lod])
  } else {
    bases_lod <- character()
    canon_lod <- character()
  }
  
  # rename to canonical forms (e.g., "IL-6_p", "IL-6.p.LOD")
  new_names <- names(DT)
  if (length(keep_val)) {
    new_names[match(keep_val, names(DT))] <- paste0(canon_val, "_", type)
  }
  if (length(keep_lod)) {
    new_names[match(keep_lod, names(DT))] <- paste0(canon_lod, ".", type, ".LOD")
  }
  setnames(DT, old = names(DT), new = new_names)
  
  # --- reshape to 4 columns ---
  # identify (renamed) value/lod columns
  val_cols <- grep(paste0("_", type, "$"), names(DT), value = TRUE)
  lod_cols <- grep(paste0("\\.", type, "\\.LOD$"), names(DT), value = TRUE)
  
  # cytokine stems (canonical bases) for alignment
  stems_val <- sub(paste0("_", type, "$"), "", val_cols)
  
  if (length(lod_cols)) {
    stems_lod <- sub(paste0("\\.", type, "\\.LOD$"), "", lod_cols)
    # align LODs to values
    lod_cols <- lod_cols[match(stems_val, stems_lod)]
    
    # safety checks
    stopifnot(length(val_cols) == 20L, length(lod_cols) == 20L)
    stopifnot(identical(stems_val, sub(paste0("\\.", type, "\\.LOD$"), "", lod_cols)))
    
    # melt paired
    long <- melt(
      DT,
      id.vars = "ppt_id",
      measure = list(value = val_cols, undetected = lod_cols),
      variable.name = "cytokine",
      value.name = c("value", "undetected")
    )
    long[, cytokine := stems_val[cytokine]]
    
  } else {
    # no LODs for this type → melt values only, fill detected with NA
    long <- melt(
      DT,
      id.vars = "ppt_id",
      measure.vars = val_cols,
      variable.name = "cytokine",
      value.name = "value"
    )
    long[, cytokine := sub(paste0("_", type, "$"), "", cytokine)]
    long[, detected := NA_real_]
  }
  
  setcolorder(long, c("ppt_id", "cytokine", "value", "undetected"))
  setorder(long, ppt_id, cytokine)
  
  # optional size checks (uncomment if you expect exactly 20 cytokines per subject)
  # stopifnot(uniqueN(long$cytokine) == 20L)
  # stopifnot(nrow(long) == 20L * uniqueN(DT$ppt_id))
  
  return(long[])
}

# "value" is raw scale, "value_new" is log transform + z-score
violin_plot <- function(dt, var = "value"){
  ggplot(dt, aes(x = "", y = .data[[var]])) +
    geom_violin(trim = FALSE, fill = "lightblue", alpha = 0.6) +
    geom_jitter(width = 0.1, alpha = 0.5, size = 0.5) +
    facet_wrap(~ cytokine, scales = ifelse(var == "value", "free_y", "fixed")) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.text.x  = element_blank(),
          axis.ticks.x = element_blank()) +
    ylab("pg/mL") + 
    theme(text = element_text(size = 9)) 
}

histogram_plot <- function(dt, var = "value"){
  ggplot(dt, aes(x = .data[[var]])) +
    geom_histogram(bins = 30, col = "black") +
    facet_wrap(~ cytokine, scales = ifelse(var == "value", "free_x", "fixed")) +
    theme_bw() +
    xlab("pg/mL") + ylab("Count") + 
    theme(text = element_text(size = 9)) 
}

# Turn alpha, gamma, and beta into Greek letters for plots.
greekify <- function(x) {
  x <- gsub("-Alpha", "\u2011\u03B1", x, fixed = TRUE)  # non-breaking hyphen + α
  x <- gsub("-Gamma", "\u2011\u03B3", x, fixed = TRUE)  # γ
  x <- gsub("IL-1Beta", "IL-1\u03B2", x, fixed = TRUE)  # β
  x
}

append_suffix <- function(dt, suffix, exclude = "ppt_id") {
  stopifnot(is.data.table(dt))
  stopifnot(is.character(suffix), length(suffix) == 1)
  
  # Find columns to rename
  cols_to_rename <- setdiff(names(dt), exclude)
  
  # Construct new names
  new_names <- paste0(cols_to_rename, "_", suffix)
  
  # Apply renaming in place
  setnames(dt, old = cols_to_rename, new = new_names)
  
  return(invisible(dt))
}

# Read in data ------------------------------------------------------------
cyto_m <- reshape_cytokines("m")
cyto_p <- reshape_cytokines("p")
cyto_c <- reshape_cytokines("c")

cyto_m[, value_new := scale(log(value)), by = cytokine]
cyto_p[, value_new := scale(log(value)), by = cytokine]
cyto_c[, value_new := scale(log(value)), by = cytokine]

cyto_m[, cytokine := greekify(cytokine)]
cyto_p[, cytokine := greekify(cytokine)]
cyto_c[, cytokine := greekify(cytokine)]

covariates <- fread(covariate_data_path)

cov <- covariates[, .(ppt_id, 
                      mat_age, 
                      bmi.procedure,
                      parity,
                      mat_educ,
                      gest.age.procedure,
                      grav,
                      marital, # 4 missing
                      mat_race_eth, # 11 missing
                      foreign, # 37 missing
                      employed, # 2 missing
                      insurance, # 1 missing
                      gest_multiple,
                      smoke, # 12 missing
                      drugs # 1 missing
                  )]

# Violin Plots (9/10 TH: Later realized these give pretty misleading ideas of the distribution) ------------------------------------------------------------
# violin_plot(cyto_p) + 
#   ggtitle("Placental cytokines (raw scale)")
# ggsave(here(results_path, "Violin Plot of Placental Cytokines - Raw Scale.png"))
# 
# violin_plot(cyto_p, "value_new") + 
#   ggtitle("Placental cytokines (log transform then z-score)")
# ggsave(here(results_path, "Violin Plot of Placental Cytokines - Log Transform and Z score.png"))
# 
# violin_plot(cyto_c) + 
#   ggtitle("Cord blood cytokines (raw scale)")
# ggsave(here(results_path, "Violin Plot of Cord Blood Cytokines - Raw Scale.png"))
# 
# violin_plot(cyto_c, "value_new") + 
#   ggtitle("Cord blood cytokines (log transform then z-score)")
# ggsave(here(results_path, "Violin Plot of Cord Blood Cytokines - Log Transform and Z score.png"))

# Histograms --------------------------------------------------------------
histogram_plot(cyto_m) + 
  ggtitle("Maternal cytokines (raw scale)")
ggsave(here(results_path, "Histogram Plot of Maternal Cytokines - Raw Scale.png"))

histogram_plot(cyto_m, "value_new") + 
  ggtitle("Maternal cytokines (log transform then z-score)")
ggsave(here(results_path, "Histogram Plot of Maternal Cytokines - Log Transform and Z score.png"))

histogram_plot(cyto_p) + 
  ggtitle("Placental cytokines (raw scale)")
ggsave(here(results_path, "Histogram Plot of Placental Cytokines - Raw Scale.png"))

histogram_plot(cyto_p, "value_new") + 
  ggtitle("Placental cytokines (log transform then z-score)")
ggsave(here(results_path, "Histogram Plot of Placental Cytokines - Log Transform and Z score.png"))

histogram_plot(cyto_c) + 
  ggtitle("Cord blood cytokines (raw scale)")
ggsave(here(results_path, "Histogram Plot of Cord Blood Cytokines - Raw Scale.png"))

histogram_plot(cyto_c, "value_new") + 
  ggtitle("Cord blood cytokines (log transform then z-score)")
ggsave(here(results_path, "Histogram Plot of Cord Blood Cytokines - Log Transform and Z score.png"))

# Data --------------------------------------------------------------------
cyto_m_wide <- dcast(cyto_m, ppt_id ~ cytokine, value.var = "value_new")
cyto_m_cor  <- cor(cyto_m_wide[, 2:21])
cyto_p_wide <- dcast(cyto_p, ppt_id ~ cytokine, value.var = "value_new")
cyto_p_cor  <- cor(cyto_p_wide[, 2:21])
cyto_c_wide <- dcast(cyto_c, ppt_id ~ cytokine, value.var = "value_new")
cyto_c_cor  <- cor(cyto_c_wide[, 2:21], use = "complete.obs")

png(here(results_path, "Correlation Matrix - Maternal.png"), width = 2000, height = 2000, res = 300)
corrplot::corrplot(cyto_m_cor, method = "color",
         addCoef.col = "black",    # print correlations on top
         number.digits = 2, number.cex = 0.6,
         tl.col = "black", tl.cex = 0.7,
         type = "upper",
         diag = FALSE)
dev.off()

png(here(results_path, "Correlation Matrix - Placental.png"), width = 2000, height = 2000, res = 300)
corrplot::corrplot(cyto_p_cor, method = "color",
         addCoef.col = "black",    # print correlations on top
         number.digits = 2, number.cex = 0.6,
         tl.col = "black", tl.cex = 0.7,
         type = "upper",
         diag = FALSE)
dev.off()

png(here(results_path, "Correlation Matrix - Cord Blood.png"), width = 2000, height = 2000, res = 300)
corrplot::corrplot(cyto_c_cor, method = "color",
         addCoef.col = "black",    # print correlations on top
         number.digits = 2, number.cex = 0.6,
         tl.col = "black", tl.cex = 0.7,
         type = "upper",
         diag = FALSE)
dev.off()

# Merging in covariates --------------------------------------------------------------

# Modify data.table wide in-place.
append_suffix(cyto_m_wide, "m")
append_suffix(cyto_p_wide, "p")
append_suffix(cyto_c_wide, "c")

# Merge back maternal, placental, and cord blood cytokines, and maternal covariates. 
dt <- cyto_m_wide |>
  merge(cyto_p_wide, by = "ppt_id") |>
  merge(cyto_c_wide, by = "ppt_id") |>
  merge(cov, by = "ppt_id")

dt[, parity := factor(parity, levels = c(0, 1), labels = c("No prior births", ">=1 prior births"))]
dt[, mat_educ := factor(mat_educ, levels = c(1,2,3), labels = c("Less than HS or HS", "Some college", "College or graduate"))]
dt[, marital := factor(marital, levels = c(1,2), labels = c("Single", "Married or living with partner"))]
dt[, foreign := factor(foreign, levels = c(0,1), labels = c("Born in US", "Born outside US"))]
dt[, gest_multiple := factor(gest_multiple, levels = c(0,1), labels = c("Singleton", "Multiples"))]
dt[, employed := factor(employed, levels = c(0, 1), labels = c("Unemployed", "employed"))]
dt[, insurance := factor(insurance, levels = c(1, 2, 3), labels = c("Private", "Public", "Self-Pay"))]
dt[, drugs := factor(drugs, level = c(0,1,2), labels = c("Never user", "Former User", "Current User"))]
dt[, smoke := factor(smoke, level = c(0,1,2), labels = c("Never smoker", "Former smoker", "Current smoker"))]

fwrite(dt, here(final_data_path, "full(n=106)_cytokine_mpc_log-z-transform_w_covariates.csv"))
