library(data.table)
library(prodlim)
library(riskRegression)
library(mstate)

source("R/data-generation.R")
source("R/crprep-timefixed.R")

set.seed(598563)

params <- list(
  "cause1" = list(
    "formula" = ~ X + Z,
    "betas" = c(0.75, 0),
    "p" = 0.25,
    "base_rate" = 1,
    "base_shape" = 0.75
  ),
  "cause2" = list(
    "formula" = ~ X + Z,
    "betas" = c(0.75, 0),
    "base_rate" = 1,
    "base_shape" = 0.75
  )
)

df_raw <- generate_dataset(
  n = 2500,
  args_event_times = list(
    mechanism = "correct_FG",
    params = params,
    censoring_type = "exponential",
    censoring_params = list("exponential" = 0.1)
  ),
  args_missingness = list(mech_params = list("prob_missing" = 0.4, "mechanism_expr" = "Z"))
)

table(df_raw$D)

# Divide cens into smaller intervals?
# Gauss quadrature? nah probs just some intervals, on quantiles of event times
# but can use outer
long_dat <- crprep.timefixed(
  Tstop = "time",
  status = "D",
  data = data.frame(df_raw),
  id = "id",
  trans = 1, # Check that this is general
  keep = c("X", "Z")
)

# Complete cases
mod <- coxph(
  Surv(Tstart, Tstop, status == 1) ~ X + Z,
  data = long_dat,
  weights = weight.cens,
  control = survival::coxph.control(timefix = FALSE)
)

basehaz(mod) |> plot()
summary(mod)

# USE SURVSPLIT!!!

# Test Beesley ------------------------------------------------------------




# Test actual smcfcs ------------------------------------------------------


# Using weighted data
