library(data.table)
library(here)
library(MAVE)
library(SuperLearner)
library(tmle)
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

# Source functions --------------------------------------------------------
source(here("R/csPCA.R"))
source(here("R/gcomp.R"))

#--- file paths ----------------------------------------------------------------
set.seed(123)
DT     <- fread(here("data/cytokines/processed/analysis_maternal-to-cord.csv"),
                stringsAsFactors = TRUE, na.strings = c("", NA))
y_name <- "IL-6_c" 
x_cols <- 2:21
c_cols <- 24:32


Y <- DT[[y_name]]
X <- DT[, x_cols, with = FALSE]
C <- DT[, c_cols, with = FALSE] 

SL.library <- c("SL.glm", 
                "SL.glmnet", 
                "tmle.SL.dbarts2")

fit <- csPCA(
  Y         = DT[[y_name]],
  X         = DT[, x_cols, with = FALSE],
  C         = DT[, c_cols, with = FALSE],
  mu_args   = list(SL.library = SL.library, method = "method.NNLS2")
)
  

# Cross-validation
mavedim <- MAVE::mave.dim( fit$mave )
cv_table <- data.frame(d=1:10, cv = mavedim$cv)
d0 <- mavedim$dim.min + 1

beta3 <- fit$mave$dir[[3]]

Z <- as.matrix(X) %*% beta3

set.seed(123)
muhat <- gcomp(Y, Z, C, SL.library = SL.library)

plot(Y, muhat)


plot(Z[,3], muhat)

newdf <- cbind(muhat_Z = muhat, 
               muhat_X = fit$mu_X, 
               Y, Z, X)


ggplot(newdf, aes(muhat_X, Y - muhat_X)) + 
  geom_point() + 
  geom_smooth()

ggplot(newdf, aes(`IL-12_m`, muhat_Z)) + 
  geom_point() + 
  geom_smooth(method = "lm")

plot(muhat, Y)
abline(0,1)

# SDR analysis ------------------------------------------------------------


sdr <- MAVE::mave(Y ~ ., data = cbind(Y, X), method = "meanMAVE")#$dir[[3]]

beta1  <- fit$mave$dir[[1]]                      
beta2  <- fit$mave$dir[[2]]                      
beta3  <- fit$mave$dir[[3]]                      
beta4  <- fit$mave$dir[[4]]                      
beta5  <- fit$mave$dir[[5]]                      
    
# Orthonormalize (preserve rownames) then varimax
Q <- qr.Q(qr(beta3))
Q <- qr.Q(qr(sdr))
rownames(Q) <- rownames(beta)
beta_final <- varimax(Q)$loadings               # "loadings" class
beta_final <- unclass(beta_final)               # plain matrix

beta.f <- data.frame(beta_final)








# Clean row names:
# 1. remove surrounding backticks if present
# 2. remove "_m" suffix
clean_names <- rownames(beta.f)
clean_names <- gsub("^`|`$", "", clean_names)   # strip leading/trailing backticks
clean_names <- sub("_m$", "", clean_names)      # drop trailing _m
rownames(beta.f) <- clean_names
# Convert to long format
df_long <- data.frame(
  Cytokine = rep(rownames(beta.f), times = ncol(beta.f)),
  Factor   = rep(colnames(beta.f),  each = nrow(beta.f)),
  Loading  = as.vector(as.matrix(beta.f))
)

# Label text = actual loadings
df_long$label <- sprintf("%.2f", df_long$Loading)

# Dynamic text color: white on dark cells, black otherwise
df_long$text_color <- ifelse(abs(df_long$Loading) > 0.4, "white", "black")

# Make Factor an ordered factor, reverse levels so X1 is on top
df_long$Factor <- factor(df_long$Factor, levels = rev(colnames(beta.f)))

# Now plot
ggplot(df_long, aes(x = Cytokine, y = Factor, fill = Loading)) +
  geom_tile(color = "grey80") +
  geom_text(aes(label = label, color = text_color), size = 3) +
  scale_color_identity() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  coord_fixed() +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_minimal(base_size = 12) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid = element_blank(),
    legend.position = "none"
  )

ggsave(filename = here("results/factor_loadings.png"), width = 11, height = 3)



# Plotting Z1, Z2, Z3 vs. the estimated muhat -----------------------------
muhat <- gcomp(Y, as.matrix(X) %*% beta_final, C, SL.library = SL.library)

