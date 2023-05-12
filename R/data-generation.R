# For simulating from Weibull distribution in shape + rate parametrization
rweibull_KM <- function(n, shape, rate) {
  (-log(runif(n)) / rate)^(1 / shape)
}

# Generate covariates
generate_covariates <- function(n, X_type) {
  dat <- data.table(id = seq_len(n), Z = rnorm(n = n, mean = 0, sd = 1))
  dat[, X := switch(
    X_type,
    "binary" = factor(rbinom(n = n, size = 1L, prob = plogis(Z))),
    "normal" = rnorm(n = n, mean = 0.5 * Z, sd = 1)
  )]
  return(dat)
}

# Generate cause 1 times in the 'squeezing' mechanism/"correct_FG"
generate_direct_times <- function(U, x, params) {
  p <- params$p
  hr <- exp(x %*% params$betas)
  val <- -log(1 - (1 - (1 - U * (1 - (1 - p)^hr))^(1 / hr)) / p)
  (val / params$base_rate)^(1 / params$base_shape)
}

# Generate event times
add_event_times <- function(dat,
                            mechanism = c("correct_FG", "misspec_FG"),
                            params,
                            censoring_params = list(
                              "exponential" = 0.5,
                              "curvy_uniform" = c(0.5, 5),
                              "curvyness" = 0.3
                            ),
                            censoring_type = c("none", "exponential", "curvy_uniform")) {

  # Match arguments
  mechanism <- match.arg(mechanism)
  censoring_type <- match.arg(censoring_type)
  predictor_formulas = list(
    "cause1" = params$cause1$formula,
    "cause2" = params$cause2$formula
  )

  # Prepare model matrices
  modmats <- lapply(predictor_formulas, function(form) {
    x <- model.matrix(form, data = dat)
    x[, !colnames(x) %in% "(Intercept)"]
  })
  x_cause1 <- modmats$cause1
  x_cause2 <- modmats$cause2
  n <- nrow(dat)

  # Generate event times
  if (mechanism == "correct_FG") {

    # Draw event indicator
    dat[, D := 1 + rbinom(
      n = .N,
      size = 1L,
      prob = (1 - params$cause1$p)^exp(x_cause1 %*% params$cause1$betas) # = P(D = 2 | X, Z)
    )]

    # Draw event times conditional on indicator
    cause1_ind <- dat[["D"]] == 1
    dat[D == 1, time := generate_direct_times(
      U = runif(.N),
      x = x_cause1[cause1_ind, ],
      params = params$cause1
    )]
    dat[D == 2, time := rweibull_KM(
      n = .N,
      shape = params$cause2$base_shape,
      rate = params$cause2$base_rate * exp(x_cause2[!cause1_ind, ] %*% params$cause2$betas)
    )]

  } else if (mechanism == "misspec_FG") {

    # Draw standard latent times
    time_cause1 <- rweibull_KM(
      n = n,
      shape = params$cause1$base_shape,
      rate = params$cause1$base_rate * exp(x_cause1 %*% params$cause1$betas)
    )

    time_cause2 <- rweibull_KM(
      n = n,
      shape = params$cause2$base_shape,
      rate = params$cause2$base_rate * exp(x_cause2 %*% params$cause2$betas)
    )

    # Pick minimum
    dat[, ':=' (
      D = 1 + as.numeric(time_cause2 < time_cause1),
      time = pmin(time_cause1, time_cause2)
    )]
  }

  # Add censoring
  if (censoring_type %in% c("exponential", "curvy_uniform")) {

    # Draw censoring times
    dat[, cens_time := switch(
      censoring_type,
      exponential = rexp(.N, rate = censoring_params$exponential),
      curvy_uniform = (censoring_params$curvy_uniform[1] - censoring_params$curvy_uniform[2]) *
        runif(.N)^(1 / censoring_params$curvyness) + censoring_params$curvy_uniform[2]
    )]

    dat[cens_time < time, ':=' (D = 0, time = cens_time)]
  }

  # Make sure D is a factor for the imputations
  dat[, D := factor(D)]
  setorder(dat, "time")

  return(dat)
}

# Add missing values
add_missingness <- function(dat,
                            mech_params = list("prob_missing" = 0.4, "mechanism_expr" = "Z")) {

  if (mech_params$prob_missing > 0) {

    # Compute part of the linear predictor conditional on covariate
    #  (mechanism_expr could also be e.g. "0.5 * Z - X" = mix of MAR and MNAR)
    linpred_expr <- parse(text = mech_params$mechanism_expr)
    linpred <- eval(linpred_expr, envir = dat)

    # Shift intercept such that average prop missing = prob_missing
    intercept_shift <- uniroot(
      f = function(eta_0, linpred, prob) mean(plogis(eta_0 + linpred)) - prob,
      interval = c(-25, 25),
      extendInt = "yes",
      linpred = linpred,
      p = mech_params$prob_missing
    )$`root`

    # Generate missing indicator, and add NAs
    dat[, ':=' (
      X_missind = rbinom(n = nrow(dat), size = 1, prob = plogis(intercept_shift + linpred)),
      X_obs = X
    )]
    dat[X_missind == 1, X := ifelse(is.numeric(X), NA_real_, NA_character_)]
  }

  return(dat)
}

# A tailored version of mice::nelsonaalen()
compute_marginal_cumhaz <- function(timevar,
                                    statusvar,
                                    cause, # number/character indicating cause of interest
                                    timefix = FALSE,
                                    type = c("cause_spec", "subdist")) {

  type <- match.arg(type)

  # First part is just like mice::nelsonaalen()
  # Second part yields cumulative subdistribution hazard
  if (type == "cause_spec") {
    mod <- survival::coxph(
      Surv(time, status) ~ 1,
      control = survival::coxph.control(timefix = timefix),
      data = cbind.data.frame("time" = timevar, "status" = as.numeric(statusvar == cause))
    )
  } else if (type == "subdist") {
    long_dat <- crprep.timefixed(
      Tstop = "time",
      status = "status",
      data = cbind.data.frame("time" = timevar, "status" = statusvar),
      trans = cause # Check that this is general
    )
    mod <- coxph(
      Surv(Tstart, Tstop, status == 1) ~ 1,
      data = long_dat,
      weights = weight.cens,
      control = survival::coxph.control(timefix = timefix)
    )
  }

  # Return cumulative hazard
  hazard <- survival::basehaz(mod)
  idx <- match(timevar, hazard[, "time"])
  return(hazard[idx, "hazard"])
}

# Find better way to use this; all code up until here was cause number agnostic
add_cumhaz_to_dat <- function(dat) {

  dat[, ':=' (
    H_cause1 = compute_marginal_cumhaz(
      timevar = time,
      statusvar = D,
      cause = 1,
      type = "cause_spec"
    ),
    H_cause2 = compute_marginal_cumhaz(
      timevar = time,
      statusvar = D,
      cause = 2,
      type = "cause_spec"
    ),
    H_subdist_cause1 = compute_marginal_cumhaz(
      timevar = time,
      statusvar = D,
      cause = 1,
      type = "subdist"
    ), # Add also cause 1 indicator
    D_star = as.numeric(D == 1)
  )]

  return(dat)
}

# Wrapper function to generate single dataset
generate_dataset <- function(n,
                             args_event_times,
                             args_missingness,
                             args_covariates = list(X_type = "binary")) {

  dat <- do.call(generate_covariates, args = c(list(n = n), args_covariates))
  do.call(add_event_times, args = c(list(dat = dat), args_event_times))
  do.call(add_missingness, args = c(list(dat = dat), args_missingness))
  return(dat[])
}

# On the large dataset, estimate the least-false parameters to be used as
# data-generating values is the misspecified FG setting
recover_weibull_lfps <- function(large_dat,
                                 params_correct_FG) {

  # Note: survreg not affected by timefix issues, see vignette

  # Recover LFPs for cause 1
  form_cs1 <- update(params_correct_FG$cause1$formula, Surv(time, D == 1) ~ .)
  mod_weib_cs1 <- survreg(form_cs1, data = large_dat, dist = "weibull")
  shape_cs1 <- 1 / mod_weib_cs1$scale
  base_rate_cs1 <- exp(mod_weib_cs1$coefficients[[1]])^(-shape_cs1)
  lfps_cs1 <- -mod_weib_cs1$coefficients[-1] * shape_cs1

  # Recover LFPs for cause 1
  form_cs2 <- update(params_correct_FG$cause2$formula, Surv(time, D == 2) ~ .)
  mod_weib_cs2 <- survreg(form_cs2, data = large_dat, dist = "weibull")
  shape_cs2 <- 1 / mod_weib_cs2$scale
  base_rate_cs2 <- exp(mod_weib_cs2$coefficients[[1]])^(-shape_cs2)
  lfps_cs2 <- -mod_weib_cs2$coefficients[-1] * shape_cs2

  # Put them in the right format
  params <- list(
    "cause1" = list(
      "formula" = params_correct_FG$cause1$formula,
      "betas" = lfps_cs1,
      "base_rate" = base_rate_cs1,
      "base_shape" = shape_cs1
    ),
    "cause2" = list(
      "formula" = params_correct_FG$cause2$formula,
      "betas" = lfps_cs2,
      "base_rate" = base_rate_cs2,
      "base_shape" = shape_cs2
    )
  )

  return(params)
}

# To determine 'true' (least-false) regression coefficients in the misspecified FG settings
recover_FG_lps <- function(large_dat,
                           censoring_type,
                           params) {

  # Predictors part of the formula
  form_rhs <- params$cause1$formula

  if (censoring_type == "none") {
    # We have to set all event 2 times to large number, larger then max event 1 time
    max_ev1_time <- large_dat[D == 1, .(time = max(time))][["time"]]
    eps <- 0.1
    large_dat[D == 2, time := max_ev1_time + eps]
    mod <- coxph(update(form_rhs, Surv(time, D == 1) ~ .), data = large_dat)
    coefs <- coef(mod)

    # If censoring, we can make use of the the fact that we know the censoring times
    # this saves looaaaads of comp time
  } else {
    large_dat[, "time_star" := ifelse(D == 2, cens_time, time)]
    mod <- coxph(update(form_rhs, Surv(time_star, D == 1) ~ .), data = large_dat)
    coefs <- coef(mod)
  }

  # Put into a data.frame, which will be easier for targets
  res <- data.frame(
    "term" = names(coefs),
    "coefs" = unname(coefs),
    "censoring_type" = censoring_type
  )
  return(res)
}
