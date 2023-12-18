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

# Compute true cumulative incidence of cause 1
compute_true_cuminc <- function(t,
                                newdat,
                                params,
                                model_type = c("correct_FG", "misspec_FG")) {

  # Make long version of newdat, and get model matrices with it
  newdat_long <- newdat[rep(seq_len(nrow(newdat)), each = length(t)), ]
  rownames(newdat_long) <- NULL

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
  cuminc <- if (model_type == "misspec_FG") {

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

    mapply(
      FUN = function(t, id, ...) integrate(prod, lower = 0, upper = t, id = id, ...)$value,
      t = t, # t here will get recycled appropriately
      id = seq_len(nrow(newdat_long)),
      MoreArgs = list(cause = 1)
    )

  } else if (model_type == "correct_FG"){

    hr_sd <- drop(exp(x_cause1 %*% params[["cause1"]][["betas"]]))
    p <- params[["cause1"]][["p"]]
    shape <- params[["cause1"]][["base_shape"]]
    rate <- params[["cause1"]][["base_rate"]]

    # Direct formula
    1 - (1 - p * (1 - exp(-rate * t^shape)))^hr_sd
  }

  cbind(newdat_long, "time" = t, "cuminc" = cuminc)
}
