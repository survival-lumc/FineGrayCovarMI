# Load important objects
tar_load(reference_patients)
newdat <- cbind.data.frame(X = reference_patients$X, Z = reference_patients$Z)
failure_time_model_dyn <- "misspec_FG" # "correct_FG"#
params <- switch(
  failure_time_model_dyn,
  "correct_FG" = tar_read(true_params_correct_FG_0.65),
  "misspec_FG" = tar_read(params_weibull_lfps_0.65)
)
model_type = failure_time_model_dyn
t <- pred_timepoints


compute_true_new(
  t = pred_timepoints,
  newdat = newdat[4, ],
  params = params,
  model_type = model_type
)

# Example with old function
plot(
  compute_true_cuminc(
    t = pred_timepoints,
    newdat = newdat[4, ],
    params = params,
    model_type = model_type
  )$cuminc,
  compute_true_new(
    t = pred_timepoints,
    newdat = newdat[4, ],
    params = params,
    model_type = model_type
  )$cuminc
)
abline(0, 1)

plot(
  compute_true_cuminc(
    t = pred_timepoints,
    newdat = newdat[4, ],
    params = params,
    model_type = model_type
  )$time,
    compute_true_cuminc(
      t = pred_timepoints,
      newdat = newdat[4, ],
      params = params,
      model_type = model_type
    )$cuminc
)

lines(
  compute_true_new(
    t = pred_timepoints,
    newdat = newdat[4, ],
    params = params,
    model_type = model_type
  )$pred_times,
  compute_true_new(
    t = pred_timepoints,
    newdat = newdat[4, ],
    params = params,
    model_type = model_type
  )$cuminc
)




# Old!! -------------------------------------------------------------------



# Support function for the hazard
weibull_hazard <- function(t, x, params, type = "hazard") {
  lp <- drop(x %*% params[["betas"]])
  rate <- params[["base_rate"]] * exp(lp)
  if (type == "hazard") {
    rate * params[["base_shape"]] * t^(params[["base_shape"]] - 1)
  } else { # type == "cumulative" for cumulative hazard
    rate * t^params[["base_shape"]]
  }
}

# Based on https://github.com/sambrilleman/simsurv/blob/master/R/simsurv.R
# For 15 points
get_gk_points <- function() {
  list(
    "gk_xi" = c(
      -0.991455371120812639207,
      -0.949107912342758524526,
      -0.86486442335976907279,
      -0.7415311855993944398639,
      -0.5860872354676911302941,
      -0.4058451513773971669066,
      -0.2077849550078984676007,
      0,
      0.2077849550078984676007,
      0.405845151377397166907,
      0.5860872354676911302941,
      0.741531185599394439864,
      0.86486442335976907279,
      0.9491079123427585245262,
      0.991455371120812639207
    ),
    "gk_wi" = c(
      0.0229353220105292249637,
      0.063092092629978553291,
      0.10479001032225018384,
      0.140653259715525918745,
      0.1690047266392679028266,
      0.1903505780647854099133,
      0.204432940075298892414,
      0.209482141084727828013,
      0.204432940075298892414,
      0.1903505780647854099133,
      0.169004726639267902827,
      0.140653259715525918745,
      0.1047900103222501838399,
      0.063092092629978553291,
      0.0229353220105292249637
    )
  )
}

integrate_t_gk <- function(t, f, ...) {
  knots <- get_gk_points()
  wi_scaled <- outer(0.5 * t, knots$gk_wi)
  xi_scaled <- outer(0.5 * t, knots$gk_xi)
  f_xi <- f(xi_scaled + 0.5 * t, ...)
  rowSums(f_xi * wi_scaled)
}

# New tings ---------------------------------------------------------------


# Make long data, and get model matrices with it
newdat_long <- cbind(
  newdat[rep(seq_len(nrow(newdat)), each = length(pred_timepoints)), ],
  pred_times = rep(pred_timepoints, times = nrow(newdat))
)

# Prepare model matrices
predictor_formulas <- lapply(params, "[[", "formula")
modmats <- lapply(predictor_formulas, function(form) {
  x <- model.matrix(form, data = newdat_long) # key
  x[, !colnames(x) %in% "(Intercept)"]
})
x_cause1 <- modmats[["cause1"]]
x_cause2 <- modmats[["cause2"]]

# Closure here
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

#cbind(newdat_long, prod(newdat_long$pred_times, cause = 1))
integrate_t_gk(
  t = newdat_long$pred_times,
  f = prod,
  cause = 1
)



# Now with outer ----------------------------------------------------------


# Prepare model matrices
predictor_formulas <- lapply(params, "[[", "formula")
modmats <- lapply(predictor_formulas, function(form) {
  x <- model.matrix(form, data = newdat)
  x[, !colnames(x) %in% "(Intercept)"]
})
x_cause1 <- modmats[["cause1"]]
x_cause2 <- modmats[["cause2"]]


# First with the correct
hr_sd <- drop(exp(x_cause1 %*% params[["cause1"]][["betas"]]))
p <- params[["cause1"]][["p"]]
shape <- params[["cause1"]][["base_shape"]]
rate <- params[["cause1"]][["base_rate"]]
cuminc <- 1 - outer(
  X = (1 - p * (1 - exp(-rate * t^shape))),
  Y = hr_sd,
  FUN = "^"
)


# Now with other tings
#newdat = newdat[1, ]

newdat

ratez <- c(0.5, 1, 1.5)
shape_global <- 0.75
funz <- function(t, shape, rate) {
  rate * shape * t^(shape - 1)
}
pred_timepoints

integrate_to_t <- Vectorize(FUN = function(t, fun, ...) {
  if (!is.numeric(t) | t < 0) stop("t should be positive.")
  ifelse(
    test = t == 0,
    yes = 0,
    no = integrate(f = fun, lower = 0L, upper = t, ...)$value
  )
}, vectorize.args = "t")

outer(
  X = ratez,
  Y = pred_timepoints,
  FUN = function(x, y) {
    integrate_to_t(fun = funz, t = y, shape = shape_global, rate = x)
    #funz(y, shape_global, x)
    #integrate_to_t(fun = funz, t = y, shape = shape_global, rate = x)
    #integrate(funz, lower = 0, upper = y, shape = shape_global)$value
  }
)

newdat

weibull_hazard <- function(t, x, params, type = "hazard") {
  lp <- drop(x %*% params[["betas"]])
  rate <- params[["base_rate"]] * exp(lp)
  if (type == "hazard") {
    #outer(t^(params[["base_shape"]] - 1), rate * params[["base_shape"]])
    outer(rate * params[["base_shape"]], t^(params[["base_shape"]] - 1))
  } else { # type == "cumulative" for cumulative hazard
    #outer(t^params[["base_shape"]], rate)
    outer(rate, t^params[["base_shape"]])
    # Change order of outer?
  }
}

ci_func <- Vectorize(function(upp, ...) {
  integrate(weibull_hazard, lower = 0, upper = upp, ...)$value
})

weibull_hazard(t = t, x = x_cause1[1:2, ], params = params[["cause1"]])

# works
ci_func(
  upp = t,
  x = x_cause1[1, ],
  params = params[["cause1"]]
)

outer(t, x_cause1[1, ], function(x, y) {
  ci_func(
    upp = x,
    x = y,
    params = params[["cause1"]]
  )
})



cause <- 1

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

test <- prod(t, cause = 1)
#outer(t, test)
test
outer(test, test)

ci_func <- Vectorize(function(upp, ...) {
  ifelse(upp == 0, 0, integrate(prod, lower = 0, upper = upp, ...)$value)
})

outer(
  X = t,
  Y = hr_sd,
  FUN = "^"
)






# Try gauss quadrature ----------------------------------------------------







weibull_hazard <- function(t, x, params, type = "hazard") {
  lp <- drop(x %*% params[["betas"]])
  rate <- params[["base_rate"]] * exp(lp)
  if (type == "hazard") {
    #outer(t^(params[["base_shape"]] - 1), rate * params[["base_shape"]])
    outer(rate * params[["base_shape"]], t^(params[["base_shape"]] - 1))
  } else { # type == "cumulative" for cumulative hazard
    #outer(t^params[["base_shape"]], rate)
    outer(rate, t^params[["base_shape"]])
    # Change order of outer?
  }
}

t <- pred_timepoints
weibull_hazard(t, x_cause1, params[["cause1"]], type = "hazard")
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

#Works for 1 time point
test <- prod(t, cause = 1)
outer(t)

outer(t, as.matrix(newdat))

knots <- get_gk_points()
wi_scaled <- outer(0.5 * t, knots$gk_wi)
xi_scaled <- outer(0.5 * t, knots$gk_xi)

f_xi <- prod(xi_scaled + 0.5 * t, cause = 1)
f_xi * wi_scaled



# New flavour -------------------------------------------------------------

pred_timepoints
newdat_long <- cbind(
  newdat[rep(seq_len(nrow(newdat)), each = length(pred_timepoints)), ],
  pred_times = rep(pred_timepoints, times = nrow(newdat))
)

weibull_hazard <- function(t, x, params, type = "hazard") {
  lp <- drop(x %*% params[["betas"]])
  rate <- params[["base_rate"]] * exp(lp)
  if (type == "hazard") {
    rate * params[["base_shape"]] * t^(params[["base_shape"]] - 1)
  } else { # type == "cumulative" for cumulative hazard
    rate * t^params[["base_shape"]]
  }
}

# Prepare model matrices
predictor_formulas <- lapply(params, "[[", "formula")
modmats <- lapply(predictor_formulas, function(form) {
  x <- model.matrix(form, data = newdat_long) # key
  x[, !colnames(x) %in% "(Intercept)"]
})
x_cause1 <- modmats[["cause1"]]
x_cause2 <- modmats[["cause2"]]

prod(newdat_long$pred_times, cause = 1)

integrate_t_gk(
  t = newdat_long$pred_times,
  f = prod,
  cause = 1
)



# With integrate ----------------------------------------------------------


