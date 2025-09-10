library(data.table)
library(here)

source(here("R/csPCA_MAVE.R"))

# Analysis 1: Maternal to Cord Blood --------------------------------------------
MC <- fread("data/cytokines/processed/analysis_maternal-to-cord.csv", 
                  stringsAsFactors = TRUE, na.strings = c("", "NA"))

MC.cor <- cor(MC[, c(2:23)])
corrplot::corrplot(MC.cor)

res      <- csPCA(MC$`TNF‑α_c`, MC[, 2:21], C = MC[,24:32])
res.mave <- mave(dt_final$`TNF‑α_c` ~ as.matrix( dt_final[, 2:21] ))

# Analysis 2: Maternal to Placenta ----------------------------------------
MP <- fread("data/cytokines/processed/analysis_maternal-to-placental.csv", 
                  stringsAsFactors = TRUE, na.strings = c("", "NA"))

MP.cor <- cor(MP[, c(2:23)])
corrplot::corrplot(MP.cor)
res      <- csPCA(MP$`TNF‑α_p`, MP[, 2:21], C = MP[,24:32])
res      <- csPCA(MP$`IL-6_p`, MP[, 2:21], C = MP[,24:32])
res.mave <- mave(dt_final$`TNF‑α_c` ~ as.matrix( dt_final[, 2:21] ))

# Analysis 3: Placental to Cord -------------------------------------------
PC <- fread("data/cytokines/processed/analysis_placental-to-cord.csv", 
                  stringsAsFactors = TRUE, na.strings = c("", "NA"))

PC.cor <- cor(PC[, c(2:23)])
corrplot::corrplot(PC.cor)
res      <- csPCA(PC$`TNF‑α_p`, MP[, 2:21], C = MP[,24:32])
res      <- csPCA(MP$`IL-6_p`, MP[, 2:21], C = MP[,24:32])
res.mave <- mave(dt_final$`TNF‑α_c` ~ as.matrix( dt_final[, 2:21] ))