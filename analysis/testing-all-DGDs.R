# Libraries
library(data.table)
library(prodlim)
library(riskRegression)
library(survival)
library(mice)
library(ggplot2)

# Source helpers
source("R/data-generation.R")

# General parameters
t_admin_cens <- 10
n <- 2500


# DGD1: Indirect ----------------------------------------------------------


# Covariates
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

# Generate latent times
dat[, T1 := generate_time_cause1(
  params_sd_cause1 = params_sd_cause1,
  params_cs_cause2 = params_cs_cause2,
  x = .SD,
  interval_upper = 10
), by = id]

dat[, T2 := rweibull_KM(
  n = n,
  shape = params_cs_cause2[["base_shape"]],
  rate = params_cs_cause2[["base_rate"]] * exp(
    params_cs_cause2[["gamma_X"]] * X + params_cs_cause2[["gamma_Z"]] * Z
  )
)]

dat[, ':=' (time = pmin(T1, T2), D = 1 + as.numeric(T2 < T1))]
dat[time >= 10, ':=' (time = 10, D = 0)]
table(dat$D)

# Check process
FGR(Hist(time, D) ~ X + Z, cause = 1, data = dat)
coxph(Surv(time, D == 2) ~ X + Z, data = dat)


# DGD2: Standard CS-based -------------------------------------------------

params_cs_cause1 <- list(
  "base_rate" = 0.25,
  "base_shape" = 0.65,
  "gamma_X" = 0.5,
  "gamma_Z" = 0
)

params_cs_cause2 <- list(
  "base_rate" = 0.25, #0.01, #0.05,
  "base_shape" = 0.5, #0.75, # Decreasing hazard
  "gamma_X" = 0.5,
  "gamma_Z" = -0.25
)

# Latent times
dat[, T1 := rweibull_KM(
  n = n,
  shape = params_cs_cause1[["base_shape"]],
  rate = params_cs_cause1[["base_rate"]] * exp(
    params_cs_cause1[["gamma_X"]] * X + params_cs_cause1[["gamma_Z"]] * Z
  )
)]

dat[, T2 := rweibull_KM(
  n = n,
  shape = params_cs_cause2[["base_shape"]],
  rate = params_cs_cause2[["base_rate"]] * exp(
    params_cs_cause2[["gamma_X"]] * X + params_cs_cause2[["gamma_Z"]] * Z
  )
)]

dat[, ':=' (time = pmin(T1, T2), D = 1 + as.numeric(T2 < T1))]
dat[time >= 10, ':=' (time = 10, D = 0)]
table(dat$D)

# Validate
coxph(Surv(time, D == 1) ~ X + Z, data = dat)
coxph(Surv(time, D == 2) ~ X + Z, data = dat)
FGR(Hist(time, D) ~ X + Z, cause = 1, data = dat)


# DGD3: Direct (FS) -------------------------------------------------------

Z <- rnorm(n, sd = 1)
X <- rbinom(n, size = 1, prob = plogis(Z))
dat <- data.table(id = seq_len(n), X, Z)
dat[, "Z" := pmin(3, pmax(Z, -3))] # Restrict range of Z

params_cuminc_cause1 <- list(
  p = 0.15,
  base_rate = 0.3,
  base_shape = 0.6,
  beta_X = 0.5,
  beta_Z = 0.25
)

params_condit_cause2 <- list(
  base_rate = 0.3,
  base_shape = 1,
  beta_star_X = -0.25,
  beta_star_Z = 0.25
)

dat[, D := 1 + rbinom(
  n = .N,
  size = 1,
  prob = (1 - params_cuminc_cause1[["p"]])^exp(
    params_cuminc_cause1[["beta_X"]] * X + params_cuminc_cause1[["beta_Z"]] * Z
  )
)]

# Generate times
dat[D == 1, time := generate_direct_times(runif(.N), .SD, params_cuminc_cause1)]
dat[D == 2, time := rweibull_KM(
  n = .N,
  shape = params_condit_cause2[["base_shape"]],
  rate = params_condit_cause2[["base_rate"]] * exp(
    params_condit_cause2[["beta_star_X"]] * X + params_condit_cause2[["beta_star_Z"]] * Z
  )
)]

dat[D == 2, time := rexp(.N, rate = params_condit_cause2[["base_rate"]] * exp(
  params_condit_cause2[["beta_star_X"]] * X + params_condit_cause2[["beta_star_Z"]] * Z
))]

# This really affects the estimate within subset
dat[time >= 10, ':=' (time = 10, D = 0)]
table(dat$D)

# Validate
#FGR(Hist(time, D) ~ X + Z, cause = 1, data = dat)
coxph(Surv(time, D == 2) ~ X + Z, data = dat, subset = (D == 2))



# DGD4: Direct (two FGs) --------------------------------------------------


params_cuminc_cause1 <- list(
  p = 0.15,
  base_rate = 0.3,
  base_shape = 0.6,
  beta_X = 0.5,
  beta_Z = 0.25
)

params_cuminc_cause2 <- list(
  p = 0.15,
  base_rate = 0.3,
  base_shape = 0.8,
  beta_X = 0.25,
  beta_Z = 0.25
)

# Get probabilities (or draw from multinomial instead?)
dat[, ':=' (
  U_ind = runif(n), # for generating indicator
  U = runif(n), # for generating time
  p1 = 1 - (1 - params_cuminc_cause1[["p"]])^(
    exp(params_cuminc_cause1[["beta_X"]] * X + params_cuminc_cause1[["beta_Z"]] * Z)
  ),
  p2 = 1 - (1 - params_cuminc_cause2[["p"]])^(
    exp(params_cuminc_cause2[["beta_X"]] * X + params_cuminc_cause2[["beta_Z"]] * Z)
  )
)]

dat[, D := fcase(
  U_ind <= p1, 1,
  p1 < U_ind & U_ind <= p1 + p2, 2,
  U_ind > p1 + p2, 0
)]

# Generate times
dat[D == 1, time := generate_direct_times(runif(.N), .SD, params_cuminc_cause1)]
dat[D == 2, time := generate_direct_times(runif(.N), .SD, params_cuminc_cause2)]

dat[time >= 10, ':=' (time = 10, D = 0)]

