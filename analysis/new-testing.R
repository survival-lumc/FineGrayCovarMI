# Libraries
library(data.table)
library(prodlim)
library(riskRegression)
library(survival)
library(mice)
library(ggplot2)

# Source helpers
source("R/data-generation.R")


# DGD1: Indirect ----------------------------------------------------------


# Covariates
n <- 2000
Z <- rnorm(n, sd = 1)
X <- rbinom(n, size = 1, prob = plogis(Z))
dat <- data.table(id = seq_len(n), X, Z)
dat[, "Z" := pmin(3, pmax(Z, -3))] # Restrict range of Z

# Set parameters
params_sd_cause1 <- list(
  "base_rate" = 0.15,
  "base_cuminc" = 0.25, # Baseline cumulative incidence for t -> infinity
  "beta_X" = 0.75,
  "beta_Z" = 0.25
)

params_cs_cause2 <- list(
  "base_rate" = 0.25, #0.01, #0.05,
  "base_shape" = 0.5, #0.75, # Decreasing hazard
  "gamma_X" = 0.5,
  "gamma_Z" = -0.25
)


# DGD2: Standard CS-based -------------------------------------------------





# DGD3: Direct (FS) -------------------------------------------------------




# DGD4: Direct (two FGs) --------------------------------------------------


