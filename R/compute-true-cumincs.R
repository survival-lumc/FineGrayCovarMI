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

