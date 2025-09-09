# q: number of columns in C (must be >= 6 here; adjust if you want fewer)
make_Theta <- function(q, p_out = 16) {
  stopifnot(q >= 6, p_out == 16)
  Theta <- matrix(0, nrow = q, ncol = p_out)
  
  # Columns 1–5 depend mainly on C1, C2
  Theta[1, 1:5] <- c(0.8, 0.6, 0.4, 0.6, 0.5)
  Theta[2, 1:5] <- c(-0.5, 0.0, -0.3, 0.0, -0.2)
  
  # Columns 6–10 depend mainly on C3, C4 (with some small cross-terms)
  Theta[3, 6:10] <- c(0.7, 0.5, 0.6, 0.4, 0.5)
  Theta[4, 6:10] <- c(0.0, -0.4, 0.0, -0.3, 0.0)
  Theta[2, 8]    <- 0.15   # small cross-loading from C2
  Theta[1, 10]   <- -0.10  # small cross-loading from C1
  
  # Columns 11–14 depend mainly on C5, C6
  Theta[5, 11:14] <- c(0.6, 0.5, 0.0, 0.4)
  Theta[6, 11:14] <- c(-0.3, 0.0, 0.4, 0.0)
  
  # Optional: tiny noise loadings to avoid perfectly sparse columns
  # (helps some optimizers, but comment out if you want exact zeros)
  # set.seed(1); Theta <- Theta + (matrix(rnorm(q*p_out,0,0.02), q, p_out) * (Theta==0))
  
  Theta
}

make_beta <- function(type = c("first", "last", "mix")) {
  type <- match.arg(type)
  
  # --- FIRST: concentrate on X1..X4, disjoint supports across columns ---
  if (type == "first") {
    beta <- matrix(0, 20, 2)
    # col 1 on X1,X2; col 2 on X3,X4 (orthogonal by disjoint support)
    beta[1:2, 1] <- 1/sqrt(2)                   # equal weights on 1,2
    beta[3:4, 2] <- c(1, -1)/sqrt(2)            # orthogonal pattern on 3,4
    return(beta)
  }
  
  # --- LAST: sparse within X5..X20, disjoint subsets for the two columns ---
  if (type == "last") {
    beta <- matrix(0, 20, 2)
    J1 <- 5:10    # 6 coords
    J2 <- 11:16   # another 6 coords, disjoint from J1
    beta[J1, 1] <- 1/sqrt(length(J1))                   # equal weights
    beta[J2, 2] <- c(1, -1, 1, -1, 1, -1)/sqrt(length(J2))
    return(beta)
  }
  
  # --- MIX: blend first block and last block, disjoint across columns ---
  if (type == "mix") {
    beta <- matrix(0, 20, 2)
    # Column 1: X1,X2 and X5..X8 (no overlap with column 2 below)
    # Give half of the energy to {1,2} and half to {5..8}.
    beta[1:2, 1] <- 0.5                           # 2*(0.5^2) = 0.5 total
    beta[5:8, 1] <- 1/(2*sqrt(2))                 # 4*((1/(2√2))^2) = 0.5
    
    # Column 2: X3,X4 and X9..X12 (disjoint from column 1’s supports)
    beta[3:4, 2]   <- c(1, -1) * 0.5              # 0.5 energy
    beta[9:12, 2]  <- c(1, -1, 1, -1) * (1/(2*sqrt(2)))  # 0.5 energy
    return(beta)
  }
}

