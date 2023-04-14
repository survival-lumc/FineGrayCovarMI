# TO BE EDITED!!
compute_true_cuminc <- function(t,
                                newdat,
                                params,
                                model_type = c("indirect", "csh_based", "direct", "both_fg"),
                                predictor_formulas) {

  # Prepare model matrices
  modmats <- lapply(predictor_formulas, function(form) {
    x <- model.matrix(form, data = newdat)
    x[, !colnames(x) %in% "(Intercept)"]
  })
  x_cause1 <- modmats[["cause1"]]
  x_cause2 <- modmats[["cause2"]]

  # Get true values
  model_type <- match.arg(model_type)

  if (model_type == "indirect") {

    # Cause 1 first
    cumhaz_cause1 <- haz_sd_cause1(t = t, x = x_cause1, params[["cause1"]], type = "cumulative")
    cuminc_cause1 <- 1 - exp(-cumhaz_cause1)

    # Cause 2
    prod <-  function(t) {
      haz <- haz_cs_cause2(t, x_cause2, params[["cause2"]], type = "hazard")
      cumhaz_cause1 <- haz_cs_cause1(t, x_cause1, params[["cause1"]], params[["cause2"]], type = "cumulative")
      cumhaz_cause2 <- haz_cs_cause2(t, x_cause2, params[["cause2"]], type = "cumulative")

      haz * exp(-cumhaz_cause1 - cumhaz_cause2)
    }

    ci_func <- Vectorize(function(upp) {
      if (upp == 0) upp <- .Machine$double.eps
      integrate(prod, lower = 0, upper = upp)$value
    })

    res <- cbind.data.frame(
      "time" = t,
      "cause1" = cuminc_cause1,
      "cause2" = ci_func(t)
    )

  } else if (model_type == "csh_based") {

    prod <-  function(t, cause) {

      haz <- if (cause == 1) {
        haz_cs_cause2(t, x_cause1, params[["cause1"]], type = "hazard") # same func
      } else {
        haz_cs_cause2(t, x_cause2, params[["cause2"]], type = "hazard")
      }

      cumhaz_cause1 <- haz_cs_cause2(t, x_cause1, params[["cause1"]], type = "cumulative")
      cumhaz_cause2 <- haz_cs_cause2(t, x_cause2, params[["cause2"]], type = "cumulative")

      haz * exp(-cumhaz_cause1 - cumhaz_cause2)
    }

    ci_func <- Vectorize(function(upp, ...) {
      if (upp == 0) upp <- .Machine$double.eps
      integrate(prod, lower = 0, upper = upp, ...)$value
    })

    res <- cbind.data.frame(
      "time" = t,
      "cause1" = ci_func(t, cause = 1),
      "cause2" = ci_func(t, cause = 2)
    )

  } else if (model_type == "direct") {

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
      "cause1" = cuminc_cause1,
      "cause2" = cuminc_cause2
    )

  } else if (model_type == "both_fg") {

    # Cause 1
    hr_sd1 <- drop(exp(x_cause1 %*% params[["cause1"]][["betas"]]))
    p1 <- params[["cause1"]][["p"]]
    shape1 <- params[["cause1"]][["base_shape"]]
    rate1 <- params[["cause1"]][["base_rate"]]
    cuminc_cause1 <- 1 - (1 - p1 * (1 - exp(-rate1 * t^shape1)))^hr_sd1

    # Cause 2
    hr_sd2 <- drop(exp(x_cause2 %*% params[["cause2"]][["betas"]]))
    p2 <- params[["cause2"]][["p"]]
    shape2 <- params[["cause2"]][["base_shape"]]
    rate2 <- params[["cause2"]][["base_rate"]]
    cuminc_cause2 <- 1 - (1 - p2 * (1 - exp(-rate2 * t^shape2)))^hr_sd2

    res <- cbind.data.frame(
      "time" = t,
      "cause1" = cuminc_cause1,
      "cause2" = cuminc_cause2
    )
  }

  return(res)
}

