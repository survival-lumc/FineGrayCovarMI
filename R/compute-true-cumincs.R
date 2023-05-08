# TO BE EDITED!!
weibull_hazard <- function(t, x, params, type = "hazard") {
  lp <- drop(x %*% params[["betas"]])
  rate <- params[["base_rate"]] * exp(lp)
  if (type == "hazard") {
    rate * params[["base_shape"]] * t^(params[["base_shape"]] - 1)
  } else { # type == "cumulative" for cumulative hazard
    rate * t^params[["base_shape"]]
  }
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


# Will return a df of true CS and subdist hazard and cumulative hazards for given covar
compute_true_hazards <- function(t,
                                 newdat,
                                 params,
                                 model_type = c("correct_FG", "misspec_FG")) {

  #browser()

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

  # Make general integration ting to integrate hazard
  integrate_hazard <- Vectorize(function(t, haz_fun, ...) {
    if (t > 0) integrate(f = haz_fun, lower = 0, upper = t, ...)$value else 0
  }, vectorize.args = "t")

  if (model_type == "misspec_FG") {

    # This part is the same as above
    prod <- function(t, cause) {

      haz <- switch(
        cause,
        "1" = weibull_hazard(t, x_cause1, params[["cause1"]], type = "hazard"),
        "2" = weibull_hazard(t, x_cause2, params[["cause2"]], type = "hazard")
      )
      cumhaz_cause1 <- weibull_hazard(t, x_cause1, params[["cause1"]], type = "cumulative")
      cumhaz_cause2 <- weibull_hazard(t, x_cause2, params[["cause2"]], type = "cumulative")
      return(haz * exp(-cumhaz_cause1 - cumhaz_cause2))
    }

    ci_func <- Vectorize(function(upp, ...) {
      ifelse(upp == 0, 0, integrate(prod, lower = 0, upper = upp, ...)$value)
    })

    # Make a function for subdist
    get_subdisthaz_misspec <- function(t, F1, F2, haz_cs1) {
      haz_cs1 * (1 - F1 - F2) / (1 - F1)
    }

    # Now collect all the good stuff
    F1 <- ci_func(t, cause = 1)
    F2 <- ci_func(t, cause = 2)
    haz_cs1 <- weibull_hazard(t, x_cause1, params[["cause1"]], type = "hazard")
    haz_cs2 <- weibull_hazard(t, x_cause2, params[["cause2"]], type = "hazard")
    cumhaz_cs1 <- weibull_hazard(t, x_cause1, params[["cause1"]], type = "cumulative")
    cumhaz_cs2 <- weibull_hazard(t, x_cause2, params[["cause2"]], type = "cumulative")
    haz_subdist1 <- get_subdisthaz_misspec(t, F1, F2, haz_cs1)
    cumhaz_subdist1 <- integrate_hazard(
      t, Vectorize(get_subdisthaz_misspec, vectorize.args = "t"),
      F1 = F1, F2 = F2,
      haz_cs1 = haz_cs1
    )

  } else if (model_type == "correct_FG") {

    # Cause1
    hr_sd <- drop(exp(x_cause1 %*% params[["cause1"]][["betas"]]))
    p <- params[["cause1"]][["p"]]
    shape <- params[["cause1"]][["base_shape"]]
    rate <- params[["cause1"]][["base_rate"]]
    F1 <- 1 - (1 - p * (1 - exp(-rate * t^shape)))^hr_sd

    # Cause 2
    px_2 <- (1 - p)^hr_sd
    hr_condit <- drop(exp(x_cause2 %*% params[["cause2"]][["betas"]]))
    shape_2 <- params[["cause2"]][["base_shape"]]
    rate_2 <- params[["cause2"]][["base_rate"]]
    cumhaz_condit <- rate_2 * hr_condit * t^shape_2
    F2 <- (1 - exp(-cumhaz_condit)) * px_2

    #
    get_subdisthaz_correct <- function(t, p, shape, rate, hr_sd) {
      nom <- p * shape * rate * exp(-rate * t^shape) * t^(shape - 1)
      denom <- 1 - p * (1 - exp(-rate * t^shape))
      hr_sd * nom / denom
    }

    haz_subdist1 <- get_subdisthaz_correct(t, p, shape, rate, hr_sd)
    cumhaz_subdist1 <- integrate_hazard(t, Vectorize(get_subdisthaz_correct, vectorize.args = "t"),
                                        p = p, shape = shape, rate = rate, hr_sd = hr_sd)
    haz_cs1 <- haz_subdist1 * (1 + F2 / (1 - F1 - F2)) # check for floating point issues?
    subdens_2 <- hr_condit * rate_2 * shape_2 * t^(shape_2 - 1) * exp(-hr_condit * rate_2 * t^shape_2)
    haz_cs2 <- subdens_2 * (1 - F1 - F2)
    # cs cumhaz after??
  }

  res <- data.table(
    "time" = t,
    haz_subdist1,
    cumhaz_subdist1,
    haz_cs1,
    haz_cs2
  )
}
