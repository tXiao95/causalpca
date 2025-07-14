# =====================================================================
#  THREE Σ[e] CASES • SINGLE γ • INTERACTION MODEL (λ, θ, α given)
# =====================================================================
library(ggplot2)
library(grid)
library(nloptr)
library(purrr)
library(dplyr)
library(ggrepel)

generate_interaction_plot <- function(beta_raw = c(1, .2),
                                      gamma_val = 5,
                                      lambda = 1,
                                      theta  = c(1, 1),
                                      alpha  = c(1, 2),
                                      mu_C   = c(1, 1),
                                      n_sim  = 5e5) {
  ## ---------- helper objects ----------------------------------------
  beta      <- beta_raw / sqrt(sum(beta_raw^2))
  theta     <- theta[1:2]                # ensure length 2
  alpha     <- alpha[1:2]
  unitvec   <- function(v) v / sqrt(sum(v^2))
  
  A      <- matrix(c(2, 1), nrow = 2)    # loads only on X2
  Sigma_C <- diag(2)                     # identity for simplicity
  
  # Σ[e] settings (row facets)
  Sigma_list <- list(
    x1_heavy = diag(c(9, 1)),
    same     = diag(c(1, 1)),
    x2_heavy = diag(c(1, 9))
  )
  
  ## ---------- simulation-based population moments -------------------
  sim_moments <- function(Sigma_e) {
    Cmat <- mvtnorm::rmvnorm(n_sim, mean = mu_C, sigma = Sigma_C)
    emat <- mvtnorm::rmvnorm(n_sim, mean = c(0, 0), sigma = Sigma_e)
    Xmat <- Cmat %*% t(A) + emat
    
    # interaction outcome
    Yvec <- Xmat %*% beta + Cmat %*% c(gamma_val, gamma_val) +
      lambda * (Cmat %*% theta) * (Xmat %*% alpha) +
      rnorm(n_sim)
    
    list(Sigma_XX = cov(Xmat),
         Sigma_XY = cov(Xmat, Yvec),
         Sigma_e  = Sigma_e)
  }
  
  ## ---------- optimiser helpers (same COBYLA trick) -----------------
  optim_unit_vec <- function(obj) {
    geq <- function(w) list("constraints" = sum(w^2) - 1,
                            "jacobian"    = 2 * w)
    w0 <- beta                             # good starting guess
    res <- nloptr(x0 = w0, eval_f = obj, eval_g_eq = geq,
                  opts = list(algorithm = "NLOPT_LN_COBYLA",
                              xtol_rel  = 1e-8, maxeval = 1000))
    unitvec(res$solution)
  }
  
  ## ---------- tidy loop over Σ[e] ----------------------------------
  plot_df <- purrr::map_dfr(names(Sigma_list), function(tag) {
    Sig_e <- Sigma_list[[tag]]
    M     <- sim_moments(Sig_e)
    Sxx   <- M$Sigma_XX; Sxy <- M$Sigma_XY
    
    ## population directions -----------------------------------------
    omega_cs   <- unitvec(beta + lambda * sum(theta * mu_C) * alpha)   # β⋆
    omega_psp  <- unitvec(Sig_e %*% omega_cs)                          # Σ_e β⋆
    omega_pls  <- unitvec(Sxy)                                         # PLS
    omega_cca  <- unitvec(solve(Sxx, Sxy))                             # CCA
    
    # doPCA objectives (scaled/unscaled) – simulated versions
    obj_scaled <- function(w)
      - (crossprod(beta, Sig_e %*% w)[1]^2) /
      (crossprod(w,    Sig_e %*% w)[1]^2)
    
    obj_unscaled <- function(w) {
      num   <- (crossprod(beta, Sig_e %*% w)[1])^2
      denom <- (crossprod(w,    Sig_e %*% w)[1])^2
      varZ  <- crossprod(w, Sxx %*% w)[1]
      -(num / denom) * varZ
    }
    
    omega_ds <- optim_unit_vec(obj_scaled)
    omega_du <- optim_unit_vec(obj_unscaled)
    
    dist <- function(w) sqrt(sum((w - beta)^2))
    
    tibble(
      SigmaCase = factor(tag, levels = c("x1_heavy", "same", "x2_heavy")),
      method = factor(c("csPCA (β⋆)",
                        "pSPCA (Σ[e] β⋆)",
                        "PLS",
                        "CCA",
                        "doPCA (scaled)",
                        "doPCA (unscaled)"),
                      levels = c("csPCA (β⋆)",
                                 "pSPCA (Σ[e] β⋆)",
                                 "PLS",
                                 "CCA",
                                 "doPCA (scaled)",
                                 "doPCA (unscaled)")),
      x = c(omega_cs[1],  omega_psp[1],  omega_pls[1],
            omega_cca[1], omega_ds[1],   omega_du[1]),
      y = c(omega_cs[2],  omega_psp[2],  omega_pls[2],
            omega_cca[2], omega_ds[2],   omega_du[2]),
      dist_to_beta = c(dist(omega_cs),  dist(omega_psp), dist(omega_pls),
                       dist(omega_cca), dist(omega_ds),  dist(omega_du))
    )
  })
  
  ## ---------- plot --------------------------------------------------
  ggplot(plot_df, aes(x = 0, y = 0)) +
    geom_segment(aes(xend = x, yend = y, colour = method), alpha = .55,
                 arrow = arrow(length = grid::unit(0.15, "cm")), linewidth = 1.05) +
    geom_text_repel(aes(x = x, y = y,
                        label = sprintf("%.2f", round(dist_to_beta, 2)),
                        colour = method),
                    size = 3, show.legend = FALSE,
                    max.overlaps = Inf, min.segment.length = 0) +
    coord_fixed(xlim = c(-1, 1), ylim = c(-1, 1)) +
    geom_abline(slope = beta_raw[2], intercept = 0, linetype = "dashed") +
    facet_wrap(~SigmaCase, labeller = labeller(
      SigmaCase = c(x1_heavy = "high x1 Σ[e]",
                    same      = "equal var Σ[e]",
                    x2_heavy  = "high x2 Σ[e]"))) +
    labs(title = sprintf("Interaction SEM • β = (%.1f, %.1f) • γ = %d",
                         beta_raw[1], beta_raw[2], gamma_val),
         subtitle = "λ = 1, θ = (1,1), α = (1,2), μ_C = (1,1)",
         x = "ω₁", y = "ω₂") +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom",
          panel.spacing   = grid::unit(1, "lines"))
}

## run three different β vectors ---------------------------------------
p_int1 <- generate_interaction_plot(beta_raw = c(1, 0.1))
p_int2 <- generate_interaction_plot(beta_raw = c(1, 1))
p_int3 <- generate_interaction_plot(beta_raw = c(1, 3))

ggsave(p_int1, filename = "figures/betaInt_highx1.png")
ggsave(p_int2, filename = "figures/betaInt_equal.png")
ggsave(p_int3, filename = "figures/betaInt_highx2.png")
