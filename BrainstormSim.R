# True model for E(mu(X) | Z)
# True model for mu(Z)

# Set the true model. 

# For every simulation experiment

# Name the experiment "NAME". Should be dated, a unique identifier, and a short name.
# Make a directory "sim-NAME". 

# Configuration File ------------------------------------------------------

# Need to be able to extend to multiple components. 

# We always assume C ~ N(0, I_q)
# We always assume \varepsilon ~ N(0,1). This is the error on Y | X, C

# n - sample size
# p - treatment size
# q - confounder size

# A - coefficient matrix
# beta - must be length 'p'. Main effect vector. 
# theta - must be length 'p'. Interaction vector
# gamma - must be length 'q'. Main effect vector. 
# alpha - must be length 'q'. Interaction vector.
# lambda - scalar. Strength of interaction.

# f - function on treatment
# g - function on confounders
# e.dist - distribution on error model for X (normal, uniform, laplace)
# confounding_type - "additive" or "multiplicative"
# interaction - value of gamma
# local - local vs. global interactions

# model - true models or estimated for csPCA and doPCA
# SL.X - SuperLearner libraries for estimating mu(X)
# SL.Z - SuperLearner libraries for estimating mu(Z)

# methods - PCA, CCA, pCCA, causal SDR, etc.