# This is a script to keep for now as quick checks for DGMs

invisible(lapply(list.files(here("R"), full.names = TRUE), source))

args_event_times <- list(
  mechanism = "correct_FG",
  #mechanism = "misspec_FG",
  #censoring_type = "exponential",
  censoring_type = "curvy_uniform",
  censoring_params = list(
    "exponential" = 0.17,
    "curvy_uniform" = c(0.5, 5),
    "curvyness" = 0.29
  ),
  #params = tar_read(params_weibull_lfps_0.15)
  params = tar_read(true_params_correct_FG_0.65)
)

# Based on 3 million obs
# 15% censored with exp rate = 0.17
# 30% censored with exp rate = 0.49
# 50% censored with exp rate = 1.44

args_missingness <- list(mech_params = list("prob_missing" = 0, "mechanism_expr" = "Z"))

dat <- generate_dataset(
  n = 3e6,
  args_event_times = args_event_times,
  args_missingness = args_missingness
)

round(prop.table(table(dat$D)), 3)



# Check betas 1 1 cumincs -------------------------------------------------


args_event_times <- list(
  mechanism = "correct_FG",
  #mechanism = "misspec_FG",
  #censoring_type = "none",
  censoring_type = "exponential",
  #censoring_type = "curvy_uniform",
  censoring_params = list(
    "exponential" = 0.25,
    "curvy_uniform" = c(0.5, 5),
    "curvyness" = 0.29
  ),
  #params = tar_read(params_weibull_lfps_0.15)
  params = list(
    "cause1" = list(
      "formula" = ~ X + Z,
      "betas" = c(1, 1),
      "p" = 0.65,
      "base_rate" = 1,
      "base_shape" = 0.75
    ),
    "cause2" = list(
      "formula" = ~ X + Z,
      "betas" = c(1, 1),
      "base_rate" = 1,
      "base_shape" = 0.75
    )
  )
)


dat <- generate_dataset(
  n = 1000000,
  args_event_times = args_event_times,
  args_missingness = args_missingness,
  args_covariates = list("X_type" = "normal")
)

prop.table(table(dat$D))

cuminc <- prodlim(Hist(time, D) ~ X, data = dat)
plot(cuminc, cause = 1, xlim = c(0, 10))
plot(cuminc, cause = 2, add = TRUE)
