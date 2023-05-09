# Weibull hazard in KM parametrization
weibull_hazard <- function(t, x, params, type = "hazard") {
  lp <- drop(x %*% params[["betas"]])
  rate <- params[["base_rate"]] * exp(lp)
  if (type == "hazard") {
    rate * params[["base_shape"]] * t^(params[["base_shape"]] - 1)
  } else { # type == "cumulative" for cumulative hazard
    rate * t^params[["base_shape"]]
  }
}

# Generic function to use for cumulative incidences and cumulative hazards
integrate_to_t <- Vectorize(FUN = function(t, fun, ...) {
  if (!is.numeric(t) | t < 0) stop("t should be positive.")
  ifelse(
    test = t == 0,
    yes = 0,
    no = integrate(f = fun, lower = 0L, upper = t, ...)$value
  )
}, vectorize.args = "t")

# Wrapper function for all
compute_true <- function(t,
                         what = c("cuminc", "hazard", "cumhazard"),
                         hazard_type = c("causespec", "subdist"),
                         cause = 1,
                         newdat,
                         params,
                         model_type = c("correct_FG", "misspec_FG")) {

  # Prepare model matrices (for now newdat is only one row), later vectorize
  predictor_formulas <- lapply(params, "[[", "formula")
  modmats <- lapply(predictor_formulas, function(form) {
    x <- model.matrix(form, data = newdat)
    x[, !colnames(x) %in% "(Intercept)"]
  })
  x_cause1 <- modmats[["cause1"]]
  x_cause2 <- modmats[["cause2"]]

  if (model_type == "misspec_FG") {

    # Create closure - we will integrate over this for cumulative incidence
    prod <- function(t, cause) {
      haz <- switch(
        cause,
        "1" = weibull_hazard(t, x_cause1, params[["cause1"]], type = "hazard"),
        "2" = weibull_hazard(t, x_cause2, params[["cause2"]], type = "hazard")
      )
      cumhaz_cause1 <- weibull_hazard(t, x_cause1, params[["cause1"]], type = "cumulative")
      cumhaz_cause2 <- weibull_hazard(t, x_cause2, params[["cause2"]], type = "cumulative")
      haz * exp(-cumhaz_cause1 - cumhaz_cause2)
    }

    # Calculate both cumulative incidences
    F1 <- integrate_to_t(fun = prod, t = t, cause = 1)
    F2 <- integrate_to_t(fun = prod, t = t, cause = 2)

    # Now we create second closure for subdistribution hazard
    get_subdisthaz_misspec <- function(t) {
      haz_cs1 <- weibull_hazard(t, x_cause1, params[["cause1"]], type = "hazard")
      F1 <- integrate_to_t(fun = prod, t = t, cause = 1)
      F2 <- integrate_to_t(fun = prod, t = t, cause = 2)
      haz_cs1 * (1 - F1 - F2) / (1 - F1)
    }

    # Calculate all hazards
    haz_cs1 <- weibull_hazard(t, x_cause1, params[["cause1"]], type = "hazard")
    haz_cs2 <- weibull_hazard(t, x_cause2, params[["cause2"]], type = "hazard")
    cumhaz_cs1 <- weibull_hazard(t, x_cause1, params[["cause1"]], type = "cumulative")
    cumhaz_cs2 <- weibull_hazard(t, x_cause2, params[["cause2"]], type = "cumulative")
    haz_subdist1 <- get_subdisthaz_misspec(t)
    cumhaz_subdist1 <- integrate_to_t(fun = get_subdisthaz_misspec, t = t)

  } else if (model_type == "correct_FG") {

    # Calculate cumulative incidences directly - here for cause 1
    hr_subdist <- drop(exp(x_cause1 %*% params[["cause1"]][["betas"]]))
    p <- params[["cause1"]][["p"]]
    shape_1 <- params[["cause1"]][["base_shape"]]
    rate_1 <- params[["cause1"]][["base_rate"]]
    F1 <- 1 - (1 - p * (1 - exp(-rate_1 * t^shape_1)))^hr_subdist

    # Cause 2
    p2_inf <- (1 - p)^hr_subdist
    hr_condit <- drop(exp(x_cause2 %*% params[["cause2"]][["betas"]]))
    shape_2 <- params[["cause2"]][["base_shape"]]
    rate_2 <- params[["cause2"]][["base_rate"]]
    cumhaz_condit <- rate_2 * hr_condit * t^shape_2
    F2 <- (1 - exp(-cumhaz_condit)) * p2_inf

    # Closure for subdistribution hazard in this mechanism
    get_subdisthaz_correct <- function(t) {
      nom <- p * shape_1 * rate_1 * exp(-rate_1 * t^shape_1) * t^(shape_1 - 1)
      denom <- 1 - p * (1 - exp(-rate_1 * t^shape_1))
      hr_subdist * nom / denom
    }

    get_cshaz_correct <- function(t, cause) {
      F1 <- 1 - (1 - p * (1 - exp(-rate_1 * t^shape_1)))^hr_subdist
      cumhaz_condit <- rate_2 * hr_condit * t^shape_2
      F2 <- (1 - exp(-cumhaz_condit)) * p2_inf
      if (cause == 1) {
        get_subdisthaz_correct(t) * (1 + F2 / (1 - F1 - F2)) # check for floating point issues?
      } else {
        subdens_2 <- hr_condit * rate_2 * shape_2 * t^(shape_2 - 1) * exp(-hr_condit * rate_2 * t^shape_2)
        subdens_2 * (1 - F1 - F2)
      }
    }

    haz_subdist1 <- get_subdisthaz_correct(t)
    cumhaz_subdist1 <- integrate_to_t(fun = get_subdisthaz_correct, t = t)
    haz_cs1 <- get_cshaz_correct(t, cause = 1)
    haz_cs2 <- get_cshaz_correct(t, cause = 2)
    cumhaz_cs1 <- integrate_to_t(fun = get_cshaz_correct, t = t, cause = 1)
    cumhaz_cs2 <- integrate_to_t(fun = get_cshaz_correct, t = t, cause = 2)
  }

  # Sort out what to return now
  data.table(
    "time" = t,
    "cuminc_1" = F1,
    "cuminc_2" = F2,
    haz_subdist1,
    haz_cs1,
    haz_cs2,
    cumhaz_subdist1,
    cumhaz_cs1,
    cumhaz_cs2
  )
}


compute_true_cuminc <- function(t,
                                newdat,
                                params,
                                model_type = c("correct_FG", "misspec_FG")) {

  predictor_formulas <- lapply(params, "[[", "formula")

  # Prepare model matrices
  modmats <- lapply(predictor_formulas, function(form) {
    x <- model.matrix(form, data = newdat)
    x[, !colnames(x) %in% "(Intercept)"]
  })
  x_cause1 <- modmats[["cause1"]]
  x_cause2 <- modmats[["cause2"]]

  # Get true values
  model_type <- match.arg(model_type)

  if (model_type == "misspec_FG") {

    prod <- function(t, cause) {

      haz <- switch(
        cause,
        "1" = weibull_hazard(t, x_cause1, params[["cause1"]], type = "hazard"),
        "2" = weibull_hazard(t, x_cause2, params[["cause2"]], type = "hazard")
      )

      cumhaz_cause1 <- weibull_hazard(t, x_cause1, params[["cause1"]], type = "cumulative")
      cumhaz_cause2 <- weibull_hazard(t, x_cause2, params[["cause2"]], type = "cumulative")

      haz * exp(-cumhaz_cause1 - cumhaz_cause2)
    }

    # If time == 0, cuminc should also be zero, so don't bother integrating
    ci_func <- Vectorize(function(upp, ...) {
      ifelse(upp == 0, 0, integrate(prod, lower = 0, upper = upp, ...)$value)
    })

    #cuminc <- try(ci_func(t, cause = 1))
    cuminc <- ci_func(t, cause = 1)
    #if (inherits(cuminc, "try-error")) {
    #  cuminc <- ci_func(t, cause = 1, rel.tol = .Machine$double.eps^.05)
    #}

    res <- cbind.data.frame(
      "time" = t,
      "cuminc" = cuminc#,
      #"cause2" = ci_func(t, cause = 2)
    )

  } else if (model_type == "correct_FG") {

    # Cause1
    hr_sd <- drop(exp(x_cause1 %*% params[["cause1"]][["betas"]]))
    p <- params[["cause1"]][["p"]]
    shape <- params[["cause1"]][["base_shape"]]
    rate <- params[["cause1"]][["base_rate"]]
    cuminc_cause1 <- 1 - (1 - p * (1 - exp(-rate * t^shape)))^hr_sd

    # Cause 2
    px_2 <- (1 - p)^hr_sd
    hr_condit <- drop(exp(x_cause2 %*% params[["cause2"]][["betas"]]))
    cumhaz_condit <- params[["cause2"]][["base_rate"]] * hr_condit * t^params[["cause2"]][["base_shape"]]
    cuminc_cause2 <- (1 - exp(-cumhaz_condit)) * px_2

    res <- cbind.data.frame(
      "time" = t,
      "cuminc" = cuminc_cause1#,
      #"cause2" = cuminc_cause2
    )

  }

  return(res)
}
