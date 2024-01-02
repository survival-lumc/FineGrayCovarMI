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

# Compute true cumulative incidence of cause 1 and 2
# along with CSH 1 and 2, SDH 1
# Partly taken from: https://github.com/survival-lumc/FineGrayDGM/blob/main/helpers.R
compute_true <- function(t,
                         newdat,
                         params,
                         model_type = c("correct_FG", "misspec_FG")) {

  # Make long version of newdat, and get model matrices with it
  newdat_long <- newdat[rep(seq_len(nrow(newdat)), each = length(t)), ]
  rownames(newdat_long) <- NULL
  times_long <- rep(t, times = nrow(newdat))

  # Prepare model matrices
  predictor_formulas <- lapply(params, "[[", "formula")
  modmats <- lapply(predictor_formulas, function(form) {
    x <- model.matrix(form, data = newdat_long)
    x[, !colnames(x) %in% "(Intercept)"]
  })
  x_cause1 <- modmats[["cause1"]]
  x_cause2 <- modmats[["cause2"]]

  model_type <- match.arg(model_type)

  # True values depend on DGM
  if (model_type == "misspec_FG") {

    prod <- function(t, cause, id) {
      haz <- switch(
        cause,
        "1" = weibull_hazard(t, x_cause1[id, ], params[["cause1"]], type = "hazard"),
        "2" = weibull_hazard(t, x_cause2[id, ], params[["cause2"]], type = "hazard")
      )
      cumhaz_cause1 <- weibull_hazard(t, x_cause1[id, ], params[["cause1"]], type = "cumulative")
      cumhaz_cause2 <- weibull_hazard(t, x_cause2[id, ], params[["cause2"]], type = "cumulative")
      haz * exp(-cumhaz_cause1 - cumhaz_cause2)
    }

    F1 <- mapply(
      FUN = function(times_long, id, ...) integrate(prod, lower = 0, upper = times_long, id = id, ...)$value,
      t = times_long,
      id = seq_len(nrow(newdat_long)),
      MoreArgs = list(cause = 1)
    )
    F2 <- mapply(
      FUN = function(times_long, id, ...) integrate(prod, lower = 0, upper = times_long, id = id, ...)$value,
      t = times_long, # t here will get recycled appropriately
      id = seq_len(nrow(newdat_long)),
      MoreArgs = list(cause = 2)
    )
    haz_cs1 <- weibull_hazard(times_long, x_cause1, params[["cause1"]], type = "hazard")
    haz_cs2 <- weibull_hazard(times_long, x_cause2, params[["cause2"]], type = "hazard")
    subdens_1 <- haz_cs1 * (1 - F1 - F2)
    subdens_2 <- haz_cs2 * (1 - F1 - F2)
    haz_subdist1 <- subdens_1 / (1 - F1)
    haz_subdist2 <- subdens_2 / (1 - F2)

  } else if (model_type == "correct_FG"){

    # Cause 1
    hr_subdist <- drop(exp(x_cause1 %*% params[["cause1"]][["betas"]]))
    p <- params[["cause1"]][["p"]]
    shape_1 <- params[["cause1"]][["base_shape"]]
    rate_1 <- params[["cause1"]][["base_rate"]]
    F1 <- 1 - (1 - p * (1 - exp(-rate_1 * times_long^shape_1)))^hr_subdist

    # Cause 2
    p2_inf <- (1 - p)^hr_subdist
    hr_condit <- drop(exp(x_cause2 %*% params[["cause2"]][["betas"]]))
    shape_2 <- params[["cause2"]][["base_shape"]]
    rate_2 <- params[["cause2"]][["base_rate"]]
    cumhaz_condit <- rate_2 * hr_condit * times_long^shape_2
    F2 <- (1 - exp(-cumhaz_condit)) * p2_inf

    # Compute both subdensities
    subdens_1 <- hr_subdist * (1 - p + p * exp(-rate_1 * times_long^shape_1))^(hr_subdist - 1) *
      p * exp(-rate_1 * times_long^shape_1) * rate_1 * shape_1 * times_long^(shape_1 - 1)

    subdens_2 <- hr_condit * rate_2 * shape_2 * times_long^(shape_2 - 1) *
      exp(-hr_condit * rate_2 * times_long^shape_2) * p2_inf # important

    # Compute all hazards
    haz_subdist1 <- subdens_1 / (1 - F1)
    haz_subdist2 <- subdens_2 / (1 - F2)
    haz_cs1 <- subdens_1 / (1 - F1 - F2)
    haz_cs2 <- subdens_2 / (1 - F1 - F2)
  }

  # Return everything in long format
  dat_ev1 <- data.table(
    newdat_long,
    "time" = times_long,
    "cause" = "1",
    "cuminc" = F1,
    "subdist_haz" = haz_subdist1,
    "cs_haz" = haz_cs1
  )
  dat_ev2 <- data.table(
    newdat_long,
    "time" = times_long,
    "cause" = "2",
    "cuminc" = F2,
    "subdist_haz" = haz_subdist2,
    "cs_haz" = haz_cs2
  )
  rbind(dat_ev1, dat_ev2)
}
