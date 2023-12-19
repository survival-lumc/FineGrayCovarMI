# Time to debug the one rep -----------------------------------------------

failure_time_model <- "misspec_FG"#"correct_FG" #
censoring_type <- "exponential" #"exponential" #  "curvy_uniform" "none"
params <- tar_read(params_weibull_lfps_0.15) #  ; both 0.65 true_params_correct_FG_0.15
tar_load(weibull_FG_lfps_0.15)
args_event_times = list(
  mechanism = failure_time_model,
  censoring_type = censoring_type,
  params = params
)
args_missingness <- list(mech_params = list("prob_missing" = prop_missing, "mechanism_expr" = "1.5 * Z"))
args_imputations <- list(m = 5, iters = 5, rjlimit = 1000) # rjlimit dont matter
args_predictions <- list(timepoints = pred_timepoints)
true_betas = switch(
  failure_time_model,
  "correct_FG" = params[["cause1"]][["betas"]],
  "misspec_FG" = weibull_FG_lfps_0.15[weibull_FG_lfps_0.15[["censoring_type"]] ==
                                            censoring_type, ][["coefs"]]
)
args_covariates = list("X_type" = "binary")
extra_args <- list(true_betas = true_betas)



# Cens expr ---------------------------------------------------------------



cens_expr <- "0.01 * exp(Z)"
censoring_type <- "exponential"
linpred_expr <- parse(text = cens_expr)
linpred <- eval(linpred_expr, envir = dat)

eval(parse(text = cens_expr), envir = dat)

dat

dat[, cens_time_test := switch(
  censoring_type,
  exponential = rexp(.N, rate = linpred),
  curvy_uniform = (censoring_params$curvy_uniform[2] - censoring_params$curvy_uniform[1]) *
    runif(.N)^(1 / censoring_params$curvyness) + censoring_params$curvy_uniform[1]
)]
