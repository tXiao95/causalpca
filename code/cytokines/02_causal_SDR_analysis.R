library(corrplot)
library(data.table)
library(here)
library(MAVE)
library(SuperLearner)

source(here("R/csPCA.R"))
source(here("R/gcomp.R"))


# Functions ---------------------------------------------------------------
correlation_matrix_plot <- function(dt, type){
  dt.cor <- cor(dt[, 2:23])
  if(type == "mc" | type == "pc"){
    highlight <- c("IL-6_c", "TNF‑α_c")
  } else{
    highlight <- c("IL-6_p", "TNF‑α_p")
  }
  # Build vectors for label colors and sizes
  label_colors <- ifelse(colnames(dt.cor) %in% highlight, "red", "black")
  label_cex    <- ifelse(colnames(dt.cor) %in% highlight, 0.7, 0.7)
  
  corrplot::corrplot(dt.cor,
                     method = "color",
                     addCoef.col = "black",    # print correlations on top
                     number.digits = 2, number.cex = 0.6,
                     tl.col = label_colors,    # customized colors
                     tl.cex = label_cex,       # customized sizes
                     type = "upper",
                     diag = TRUE, cl.pos = "n")
}

#--- robust single-analysis helper --------------------------------------------
compute_beta_final_safe <- function(DT, y_name, x_cols = 2:21, c_cols = 24:32,
                                    mu_args = list(), mave_args = list()) {
  tryCatch({
    fit <- csPCA(
      Y         = DT[[y_name]],
      X         = DT[, x_cols, with = FALSE],
      C         = DT[, c_cols, with = FALSE],
      mu_args   = mu_args,
      mave_args = mave_args
    )
    
    d     <- MAVE::mave.dim(fit$mave)$dim.min
    d <- ifelse(d==1, 2, d)
    beta  <- fit$mave$dir[[d]]                      # rows = variables, cols = directions
    
    # Orthonormalize (preserve rownames) then varimax
    Q <- qr.Q(qr(beta))
    rownames(Q) <- rownames(beta)
    beta_final <- varimax(Q)$loadings               # "loadings" class
    beta_final <- unclass(beta_final)               # plain matrix
    
    list(ok = TRUE,
         beta_final = beta_final,
         d = d,
         n_vars = nrow(beta),
         outcome = y_name,
         error = NULL)
  }, error = function(e) {
    # Return a structured failure without stopping the whole script
    list(ok = FALSE,
         beta_final = NULL,
         d = NA_integer_,
         n_vars = NA_integer_,
         outcome = y_name,
         error = conditionMessage(e))
  })
}

prelim <- function(type = c("mc", "mp", "pc")){
  
  filename <- switch (type,
    mc = "maternal-to-cord",
    mp = "maternal-to-placental",
    pc = "placental-to-cord"
  )
  
  plot_title <- switch (type,
    mc = "Correlation Matrix (with outcome) - Maternal to Cord Blood",
    mp = "Correlation Matrix (with outcome) - Maternal to Placental",
    pc = "Correlation Matrix (with outcome) - Placental to Cord Blood"
  )
  
  dt <- fread(here("data/cytokines/processed/", paste0("analysis_", filename, ".csv")), 
                    stringsAsFactors = TRUE, na.strings = c("", "NA"))
  
  png(filename = here("outputs/cytokines/figures/", paste0(plot_title, ".png")), 
      width = 2000, height = 2000, res = 200)
  correlation_matrix_plot(dt, type = type)
  dev.off()
}


# Preliminary analyses - correlation matrix -------------------------------

prelim("mc")
prelim("mp")
prelim("pc")

#--- file paths ----------------------------------------------------------------
paths <- list(
  MC = "data/cytokines/processed/analysis_maternal-to-cord.csv",
  MP = "data/cytokines/processed/analysis_maternal-to-placental.csv",
  PC = "data/cytokines/processed/analysis_placental-to-cord.csv"
)

# Outcomes per dataset (note _c vs _p)
outcomes <- list(
  MC = c("TNF‑α_c", "IL-6_c"),
  MP = c("TNF‑α_p", "IL-6_p"),
  PC = c("TNF‑α_c", "IL-6_c")
)

#--- run all 6 analyses --------------------------------------------------------
beta_final_list <- list()

for (tag in names(paths)) {
  DT <- fread(paths[[tag]], stringsAsFactors = TRUE, na.strings = c("", "NA"))
  for (y in outcomes[[tag]]) {
    key <- paste(tag, y, sep = "__")
    beta_final_list[[key]] <- compute_beta_final_safe(DT, y_name = y)
  }
}
# ensure names are set (they already are via keys)
names(beta_final_list) <- names(beta_final_list)
saveRDS(beta_final_list, here("outputs/cytokines/beta_varimax_loadings.rds"))
