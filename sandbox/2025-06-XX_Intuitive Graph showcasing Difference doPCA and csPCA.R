# 6/XX/25: This script is important! Will be use din paper
# Was used to visualize intuitively how doPCA differs from csPCA 
# using the linear model. 

# =============================================================
#  2×3 GRID — now including doPCA (scaled & unscaled)
# =============================================================
library(ggplot2)
library(grid) 
library(nloptr)
library(purrr)
library(dplyr)
library(MASS)
library(ggcorrplot)
library(tidyr)

# ---------- fixed pieces (p = 2) ------------------------------
generate_plot <- function(beta_raw = c(1, .2)){
  beta     <- beta_raw / sqrt(sum(beta_raw^2))
  
  A       <- matrix(c(2, 1), nrow = 2)      
  VarC    <- 1
  unitvec <- function(v) v / sqrt(sum(v^2))
  
  # ---------- Σ_e settings --------------------------------------
  Sigma_list <- list(
    x1_heavy   = diag(c(9,1)),
    same = diag(c(1, 1)),      # big heterosked.
    x2_heavy   = diag(c(1, 9))   # small heterosked.
  )
  
  # Sigma_XX <- A %*% VarC %*% t(A) + Sigma_e
  # ---------- γ settings ----------------------------------------
  gamma_vals <- c(5)
  
  # ---------- helper: doPCA optimiser ---------------------------
  # -------------------------------------------------------------
  # install.packages("nloptr")  # once, if not already installed
  
  # -------------------------------------------------------------
  #  Constrained optimiser for doPCA objectives (p = 2 here)
  #  * scaled  :  f_s(ω) =  (β'Σ_e ω / ω'Σ_e ω)^2
  #  * unscaled:  f_u(ω) =  (β'Σ_e ω / ω'Σ_e ω)^2 · (ω'Σ_XX ω)
  #
  #  Returns the maximiser on the unit circle.
  # -------------------------------------------------------------
  library(nloptr)
  
  optim_dopca <- function(beta, Sigma_e, Sigma_XX, scaled = TRUE) {
    ## ---- analytic gradient of the objective -------------------
    fgrad <- function(w) {
      a <- as.numeric(beta %*% Sigma_e %*% w)          # β'Σe ω   (scalar)
      b <- as.numeric(w %*% Sigma_e %*% w)          # ω'Σe ω   (scalar)
      
      if (scaled) {
        # t0 <- Sigma_e %*% beta
        # t1 <- Sigma_e %*% omega
        # t2 <- t(omega) %*% t1
        # t3 <- t(omega) %*% t0
        # t4 <- t2^4
        
        # g <- (2*t3/t2^2) * t0 - ( )
        
        
        f  <- (a^2)/(b^2)
      } else {
        wXw <- w %*% Sigma_XX %*% w         # ω'Σxx ω
        f1 <- (a^2)/(b^2)                                      # inner ratio
        f  <- f1 * wXw
      }
      return(-f)
    }
    
    ## ---- equality constraint  g(ω) = ‖ω‖² - 1 = 0 ------------
    geq <- function(w) sum(w^2) - 1
    
    ## ---- initial guess: unit vector along β ------------------
    w0 <- beta / sqrt(sum(beta^2))
    
    sol <- nloptr(
      x0        = w0,
      eval_f    = fgrad,        # returns both objective and gradient
      eval_g_eq = geq,
      opts      = list(
        algorithm = "NLOPT_LN_COBYLA",
        xtol_rel  = 1e-10,
        maxeval   = 1000
      )
    )
    
    # return(sol)
    sol$solution / sqrt(sum(sol$solution^2))     # ensure exact unit norm
  }
  
  optim_cspca <- function(beta, Sigma_XX) {
    ## ---- analytic gradient of the objective -------------------
    fgrad <- function(w) {
      a <- as.numeric(beta %*% Sigma_XX %*% w)          # β'Σe ω   (scalar)
      b <- as.numeric(w %*% Sigma_XX %*% w)          # ω'Σe ω   (scalar)
        
      f  <- t(beta) %*% Sigma_XX %*% beta - a^2 / b
      return(f)
    }
    
    ## ---- equality constraint  g(ω) = ‖ω‖² - 1 = 0 ------------
    geq <- function(w) sum(w^2) - 1
    
    ## ---- initial guess: unit vector along β ------------------
    w0 <- beta / sqrt(sum(beta^2))
    
    sol <- nloptr(
      x0        = w0,
      eval_f    = fgrad,        # returns both objective and gradient
      eval_g_eq = geq,
      opts      = list(
        algorithm = "NLOPT_LN_COBYLA",
        xtol_rel  = 1e-10,
        maxeval   = 1000
      )
    )
    
    # return(sol)
    sol$solution / sqrt(sum(sol$solution^2))     # ensure exact unit norm
  }
  
  
  # optim_dopca(beta, Sigma_e, Sigma_XX)
  # optim_dopca(beta, Sigma_e, Sigma_XX, scaled = FALSE)
  
  # ---------- tidy data for all facets --------------------------
  beta_unit <- unitvec(beta_raw)     # same function used earlier
  
  plot_df <- purrr::map_dfr(names(Sigma_list), function(sname) {
    Sigma_e <- Sigma_list[[sname]]
    
    purrr::map_dfr(gamma_vals, function(gam) {
      ## population covariances
      Sigma_XX <- A %*% VarC %*% t(A) + Sigma_e
      Sigma_XY <- Sigma_XX %*% beta_unit + A * VarC * gam
      
      ## first‐direction vectors
      omega_pca   <- unitvec( eigen(Sigma_XX)$vectors[, 1] )       # PCA
      
      if(omega_pca[1] < 0) omega_pca <- omega_pca * -1
      
      omega_cs    <- optim_cspca(beta_unit, Sigma_XX)              # pCCA & csPCA
      omega_psp   <- unitvec(Sigma_e %*% beta_unit)                # pSPCA
      omega_pls   <- unitvec(Sigma_XY)                             # PLS
      omega_cca   <- unitvec(solve(Sigma_XX, Sigma_XY))            # CCA
      omega_dopca_scaled   <- optim_dopca(beta_unit, Sigma_e, Sigma_XX, scaled = TRUE)
      omega_dopca_unscaled <- optim_dopca(beta_unit, Sigma_e, Sigma_XX, scaled = FALSE)
      
      # print(omega_dopca_unscaled)
      # print(omega_dopca_scaled)
      
      ## helper to compute Euclidean distance to β
      dist <- function(w) sqrt(sum((w - beta_unit)^2))
      
      tibble(
        SigmaCase = factor(sname, levels = c("x1_heavy", "same", "x2_heavy")),
        gamma_lab = factor(paste0("γ = ", gam),
                           levels = paste0("γ = ", gamma_vals)),
        method = factor(c("PCA", 
                          "csPCA (β)",
                          "pSPCA (Σ[e] β)",
                          "SPCA",
                          "CCA",
                          "doPCA (scaled)",
                          "doPCA (unscaled)"),
                        levels = c("PCA", 
                                   "csPCA (β)",
                                   "pSPCA (Σ[e] β)",
                                   "SPCA",
                                   "CCA",
                                   "doPCA (scaled)",
                                   "doPCA (unscaled)")),
        x = c(omega_pca[1], omega_cs[1],  omega_psp[1], omega_pls[1],
              omega_cca[1], omega_dopca_scaled[1], omega_dopca_unscaled[1]),
        y = c(omega_pca[2], omega_cs[2],  omega_psp[2], omega_pls[2],
              omega_cca[2], omega_dopca_scaled[2], omega_dopca_unscaled[2]),
        dist_to_beta = c(dist(omega_pca),
                         dist(omega_cs),
                         dist(omega_psp),
                         dist(omega_pls),
                         dist(omega_cca),
                         dist(omega_dopca_scaled),
                         dist(omega_dopca_unscaled))
      )
    })
  })
  
  
  # ---------- plot ------------------------------------------------------
  # install.packages("ggrepel")   # once, if needed
  library(ggrepel)
  
  ggplot(plot_df, aes(x = 0, y = 0)) +
    geom_segment(aes(xend = x, yend = y, colour = method), alpha = 0.5,
                 arrow = arrow(length = grid::unit(0.16, "cm")),
                 linewidth = 1.1) +
    ## distance label (rounded to 2 d.p.) --------------------------
  geom_text_repel(aes(x = x, y = y,
                      label = sprintf("%.2f", round(dist_to_beta, 2)),
                      colour = method),
                  size = 3,                  # text size
                  show.legend = FALSE,       # keep legend clean
                  min.segment.length = 0) +  # draw small leaders if needed
    coord_fixed(xlim = c(-1, 1), ylim = c(-1, 1)) +
    geom_abline(slope = beta_raw[2], intercept = 0, linetype = "dashed") +
    geom_hline(yintercept = 0, colour = "grey80") +
    geom_vline(xintercept = 0, colour = "grey80") +
    facet_wrap(~SigmaCase, labeller = labeller(
      SigmaCase = c(x1_heavy = "high x1 Σ[e]",
                    same      = "equal var Σ[e]",
                    x2_heavy = "high x2 Σ[e]"))) + 
    # facet_grid(SigmaCase ~ gamma_lab,
    #            labeller = labeller(
    #              SigmaCase = c(x1_heavy = "high x1 Σ[e]",
    #                            same      = "equal var Σ[e]",
    #                            x2_heavy = "high x2 Σ[e]"))) +
    labs(title = sprintf("Population first directions (ω) for β = (%.1f, %.1f)",
                         beta_raw[1], beta_raw[2]),
         subtitle = "Columns: Σ[e] heteroskedasticity with confounding γ=5 and A=(2,1)",
         x = "w1", y = "w2") +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom",
          panel.spacing   = grid::unit(1, "lines"))
}

p1 <- generate_plot(c(1, .1))
p3 <- generate_plot(c(1, 1))
p5 <- generate_plot(c(1, 3))

ggsave(p1, filename = "beta_high_x1.png")
ggsave(p3, filename = "beta_equal.png")
ggsave(p5, filename = "beta_high_x2.png")


# Plotting z against E(Y(z)) ----------------------------------------------

# True  beta    : 0.9950372 0.0995037
# doPCA unscaled: 0.6548505 0.7557585
# doPCA scaled:   0.3418811 0.9397432
# 
# doPCA unscaled: 0.4681053 0.8836727
# doPCA scaled  : 0.2544832 0.9670772
#
#

omega_do_unscaled <- c(0.6548505, 0.7557585)
omega_do_scaled   <- c(0.3418811, 0.9397432)

unitvec <- function(v) v / sqrt(sum(v^2))
beta_raw <- c(1, 1)
beta     <- unitvec(beta_raw)

A     <- matrix(c(2, 1), nrow = 2)      
VarC  <- 1
muC   <- 1
gamma <- 5

# ---------- Σ_e settings --------------------------------------
set.seed(12340)
Sigma_e <- diag(c(9,1))

N <- 100

e <- mvrnorm(n = N, c(0,0), Sigma_e)
C <- rnorm(n = N, muC, VarC)

X             <- C %*% t(A) + e
Z_beta        <-  X %*% beta 
Z_do_unscaled <- X %*% omega_do_unscaled
Z_do_scaled   <- X %*% omega_do_scaled

# Form data.frame with each Z against E(Y(z)) and E(Y(x))

# Causal mean E(Y(x))
mu_x <- function(x){
  return( as.numeric( beta %*% x ) )
}

mu_z <- function(x, omega = beta){
  z     <- sum(x * omega)
  num   <- t(beta) %*% Sigma_e %*% omega
  denom <- t(omega) %*% Sigma_e %*% omega
  
  return( as.numeric(  num / denom * z   ) )
}

# mu(z) vs. mu(x) plot ----------------------------------------------------

dt <- data.frame(X = X, 
                 Z_beta = Z_beta, 
                 Z_do_unscaled = Z_do_unscaled,
                 Z_do_scaled = Z_do_scaled,
                 mu_z_u = apply(X, 1, mu_z, omega = omega_do_unscaled),
                 mu_z_s = apply(X, 1, mu_z, omega = omega_do_scaled),
                 mu_x = apply(X, 1, mu_x))


mu_dat  <- dt[ , c("mu_x", "mu_z_s", "mu_z_u")]
mu_corr <- cor(mu_dat, use = "pairwise.complete.obs")

# Names you want to see in the plot
mu_labels <- c(
  "mu[x]",            # μₓ
  "mu[z]^{scaled}",   # μ_z^{scaled}
  "mu[z]^{unscaled}"  # μ_z^{unscaled}
)

GGally::ggpairs(
  mu_dat,                         # your 3-column data frame
  columnLabels = mu_labels,       # pretty strip / axis labels
  labeller      = "label_parsed", # parse as plotmath
  upper         = list(continuous = "cor"),   # (optional) show r in upper
  lower         = list(continuous = wrap(ggally_points, alpha = .7)),
  diag          = list(continuous = "densityDiag")
) +
  theme_bw()
ggsave(file = "mu(x) vs. mu(z).png")

# mu(z) vs. z plot --------------------------------------------------------

# Plot all combinations of the three Z columns on the x-axis and three mu columsn on the y-axis. 
# Reshape into long form with every Z–μ pairing
plot_dat <- dt %>% 
  pivot_longer(c(Z_beta, Z_do_unscaled, Z_do_scaled),
               names_to  = "Z_name",  values_to = "Z_val") %>% 
  pivot_longer(c(mu_x,  mu_z_s, mu_z_u),
               names_to  = "mu_name", values_to = "mu_val")

# Labels for facets (parsed so Greek letters render)
lab_Z  <- c(
  Z_beta        = "Z[beta]",
  Z_do_unscaled = "Z[do]^{unscaled}",
  Z_do_scaled   = "Z[do]^{scaled}"
)
lab_mu <- c(
  mu_x   = "mu[x]",
  mu_z_s = "mu[z]^{scaled}",
  mu_z_u = "mu[z]^{unscaled}"
)

ggplot(plot_dat, aes(x = Z_val, y = mu_val)) +
  geom_point(alpha = .6) +
  facet_grid(
    rows = vars(mu_name), cols = vars(Z_name),
    labeller = labeller(
      Z_name  = as_labeller(lab_Z,  label_parsed),
      mu_name = as_labeller(lab_mu, label_parsed)
    )
  ) +
  theme_bw() +
  labs(x = NULL, y = NULL) + 
  geom_abline(slope = 1, intercept = 0) + 
  xlab("Z (Principal Score)") + 
  ylab("Causal Mean")

ggsave(filename = "mu(z) vs z.png")
