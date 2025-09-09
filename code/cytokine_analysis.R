library(data.table)
library(forcats)
library(here)
library(naniar)

source("R/csPCA_MAVE.R")

# data --------------------------------------------------------------------

dt <- fread(here("data/Eick/full_cytokine_p_data_cleaned.csv"), 
            stringsAsFactors = TRUE, na.strings = c("", "NA"))

dt[, smoke := fct_collapse(
  smoke,
  "Former or Current smoker" = c("Current smoker", "Former smoker")
)]

dt[, drugs := fct_collapse(
  drugs,
  "Former or Current user" = c("Current User", "Former User")
)]


vis_miss(dt)
ggplot2::ggsave("results/real_data_analysis/cytokines/Missing Data Pattern.png")

drop_cols <- c("marital", "mat_race_eth", "foreign", "gest_multiple", 
               "smoke")

dt_final <- dt[, !..drop_cols]
dt_final <- dt_final[complete.cases(dt_final)] 

fwrite(dt_final, here("data/Eick/cytokine_p_data.csv"))

# csPCA -------------------------------------------------------------------
res <- csPCA(dt_final$`IL-6_c`, dt_final[, 2:21], C = dt_final[,24:32])
res.mave <- mave(dt_final$`IL-6_c` ~ as.matrix( dt_final[, 2:21] ))

a <- res$dir[[2]]
a.mave <- res.mave$dir[[4]]

vm      <- varimax(a)
vm.mave <- varimax(a.mave)
library(pheatmap)

pheatmap(vm$loadings)
         

library(ggplot2)
library(reshape2)

df_long <- melt(vm$loadings)
df_long <- melt(vm$loadingsames = c("Cytokine", "Component"), value.name = "Loading")

ggplot(df_long, aes(x = Component, y = Cytokine, fill = Loading)) +
  geom_tile(color = "grey80") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
