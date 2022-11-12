##**********************************##
## Data-generating function helpers ##
##**********************************##

# (Wrap later into generate-dataset)


# General  ----------------------------------------------------------------


#' @title Sample from Weibull distirbution in K&M parametrisation.
#'
#' @param shape Shape parameter of Weibull distribution.
#' @param rate Rate parameter of Weibull distirbution.
#' @param n Sample size
#'
#' @return n samples from wWibull distribution.
#'
rweibull_KM <- function(n, shape, rate) {
  (-log(runif(n)) / rate)^(1 / shape)
}

# Make this into more general function
add_marginal_cumhaz <- function(timevar,
                                statusvar, # competing event indicator
                                cause,
                                timefix = FALSE,
                                type = c("cause_spec", "subdist")) {

  type <- match.arg(type)

  # First part is just like mice::nelsonaalen
  if (type == "cause_spec") {
    mod <- survival::coxph(
      Surv(time, status) ~ 1,
      control = survival::coxph.control(timefix = timefix),
      data = cbind.data.frame("time" = timevar, "status" = as.numeric(statusvar == cause))
    )
  } else {
    long_dat <- crprep(
      Tstop = "time",
      status = "status",
      data = cbind.data.frame("time" = timevar, "status" = statusvar),
      trans = cause
    )
    mod <- coxph(
      Surv(Tstart, Tstop, status == 1) ~ 1,
      data = long_dat,
      weights = weight.cens,
      control = survival::coxph.control(timefix = timefix)
    )
  }

  # Get cumhaz
  hazard <- survival::basehaz(mod)
  idx <- match(timevar, hazard[, "time"])
  return(hazard[idx, "hazard"])
}


# Helpers for direct methods ----------------------------------------------


generate_direct_times <- function(U, x, params) {
  p <- params[["p"]]
  hr <- exp(x %*% params[["betas"]])
  val <- -log(1 - (1 - (1 - U * (1 - (1 - p)^hr))^(1 / hr)) / p)
  (val / params[["base_rate"]])^(1 / params[["base_shape"]])
}


# Helpers for CS/indirect methods -----------------------------------------


haz_sd_cause1 <- function(t, x, params, type = "hazard") {
  base <- log(1 / (1 - params[["base_cuminc"]]))
  lp <- drop(x %*% params[["betas"]])

  if (type == "hazard") {
    params[["base_rate"]] * exp(-params[["base_rate"]] * t / base + lp)
  } else { # Cumulative
    exp(lp) * base * (1 - exp(-params[["base_rate"]] * t / base))
  }
}

haz_cs_cause2 <- function(t, x, params, type = "hazard") {
  lp <- drop(x %*% params[["betas"]])
  rate <- params[["base_rate"]] * exp(lp)
  if (type == "hazard") {
    rate * params[["base_shape"]] * t^(params[["base_shape"]] - 1)
  } else { # Cumulative
    rate * t^params[["base_shape"]]
  }
}

# The integral to be numerically solved, see relative sd1 and cs2
integral_haz_cause1 <- Vectorize(function(t, x, params_sd_cause1, params_cs_cause2) {

  integral_fun <- function(t, x, params_sd_cause1, params_cs_cause2) {
    exponent <- haz_cs_cause2(t, x, params_cs_cause2, type = "cumulative") -
      haz_sd_cause1(t, x, params_sd_cause1, type = "cumulative")
    haz_sd_cause1(t, x, params_sd_cause1) * exp(exponent)
  }

  args <- list(
    "f" = integral_fun,
    "lower" = 0,
    "upper" = t,
    "x" = x,
    "params_sd_cause1" = params_sd_cause1,
    "params_cs_cause2" = params_cs_cause2
  )

  # Try integrating
  val <- try(expr = do.call(integrate, args), silent = TRUE)

  # t probably very large, take Inf as upper instead to avoid errors
  if (inherits(val, "try-error")) {
    args$upper <- Inf
    val_infty <- try(expr = do.call(integrate, args), silent = TRUE)

    if (inherits(val_infty, "try-error")) {
      #covars <- paste(round(c(x[["X"]], x[["Z"]]), 2), collapse = ", ")
      stop("Likely integral itself becomes infinite")
    } else val_infty$value

  } else val$value

}, vectorize.args = "t")


# Define cause-specific hazards for cause 1
haz_cs_cause1 <- function(t, x, params_sd_cause1, params_cs_cause2, type = "hazard") {
  exponent <- haz_cs_cause2(t, x, params_cs_cause2, type = "cumulative") -
    haz_sd_cause1(t, x, params_sd_cause1, type = "cumulative")
  num <- haz_sd_cause1(t, x, params_sd_cause1) * exp(exponent)

  if (any(is.infinite(num))) {
    mssg <- paste0(
      "Function lambda1(t) * exp(H2(t) - Lambda1(t)) seems to go to infinity: ",
      "try other parameter values for H2(t) such that h1(t) > 0 for all t!"
    )
    stop(mssg)
  } else if (any(is.nan(num))) {
    # R will return NaN is the number is too close too zero
    num[is.nan(num)] <- .Machine$double.eps
  }

  # Integration checks
  denom <- 1 - integral_haz_cause1(t, x, params_sd_cause1, params_cs_cause2)
  if (any(denom <= 0)) {
    lim_haz <- min(t[which(denom <= 0)])
    stop(paste0("h1(t) is undefined or negative starting at approx t = ", lim_haz))
  }

  # Return (cumulative) hazard
  if (type == "hazard") {
    num / denom
  } else {
    -log(denom)
  }
}

invert_surv_cause1 <- function(t, x, params_sd_cause1, params_cs_cause2, U) {
  H1 <- haz_cs_cause1(t, x, params_sd_cause1, params_cs_cause2, type = "cumulative")
  exp(-H1) - U
}

generate_indirect_cause1 <- function(params_sd_cause1,
                                     params_cs_cause2,
                                     x,
                                     interval_upper = 1e3) {
  # Generate u
  U <- runif(1)

  # Evaluate survival at limit
  surv_upper <- invert_surv_cause1(
    interval_upper,
    params_sd_cause1 = params_sd_cause1,
    params_cs_cause2 = params_cs_cause2,
    x = x,
    U = U
  )

  if (is.nan(surv_upper)) {
    stop("Woops!")
  } else if (surv_upper > 0) {
    interval_upper # and then censor by maxt
  } else {
    # Generate time by inversion
    uniroot(
      f = invert_surv_cause1,
      interval = c(.Machine$double.eps, interval_upper),
      extendInt = "yes",
      params_sd_cause1 = params_sd_cause1,
      params_cs_cause2 = params_cs_cause2,
      x = x,
      U = U
    )$root
  }
}


# Full wrapper ------------------------------------------------------------


generate_complete_dataset <- function(n = 2000,
                                      params,
                                      model_type = c("indirect", "csh_based", "direct", "both_fg"),
                                      predictor_formulas,
                                      X_type = c("binary", "normal"),
                                      control = list("admin_cens_time" = NULL, "cens_rate" = NULL)) {

  # Covariate generation (could even feed this as argument to other function)
  dat <- data.table(id = seq_len(n), Z = rnorm(n, mean = 0, sd = 1))

  if (X_type == "binary") {
    dat[, "X" := factor(rbinom(.N, size = 1, prob = plogis(Z)))]
  } else {
    dat[, "X" := rnorm(.N, mean = 0.5 * Z, sd = 1)]
    #dat[, "X" := pmin(3, pmax(X, -3))] # restrict range
  }

  #dat[, "Z" := pmin(3, pmax(Z, -3))] # restrict range

  # Prepare model matrices
  modmats <- lapply(predictor_formulas, function(form) {
    x <- model.matrix(form, data = dat)
    x[, !colnames(x) %in% "(Intercept)"]
  })
  x_cause1 <- modmats[["cause1"]]
  x_cause2 <- modmats[["cause2"]]

  # Generate times
  model_type <- match.arg(model_type)
  if (model_type == "indirect") {

    # Generate latent times
    dat[, T1 := generate_indirect_cause1(
      params_sd_cause1 = params[["cause1"]],
      params_cs_cause2 = params[["cause2"]],
      x = x_cause1[.I, ],
      interval_upper = control[["admin_cens_time"]]
    ), by = id]

    dat[, T2 := rweibull_KM(
      n = n,
      shape = params[["cause2"]][["base_shape"]],
      rate = params[["cause2"]][["base_rate"]] * exp(x_cause2 %*% params[["cause2"]][["betas"]])
    )]

    # Pick minimum
    dat[, ':=' (D = 1 + as.numeric(T2 < T1), time = pmin(T1, T2))]
    dat[, c("T1", "T2") := NULL]

  } else if (model_type == "csh_based") {

    # Generate latent times
    T1 <- rweibull_KM(
      n = n,
      shape = params[["cause1"]][["base_shape"]],
      rate = params[["cause1"]][["base_rate"]] * exp(x_cause1 %*% params[["cause1"]][["betas"]])
    )

    T2 <- rweibull_KM(
      n = n,
      shape = params[["cause2"]][["base_shape"]],
      rate = params[["cause2"]][["base_rate"]] * exp(x_cause2 %*% params[["cause2"]][["betas"]])
    )

    # Pick minimum
    dat[, ':=' (D = 1 + as.numeric(T2 < T1), time = pmin(T1, T2))]

  } else if (model_type == "direct") {

    # Draw event indicator
    dat[, D := 1 + rbinom(
      n = .N,
      size = 1,
      prob = (1 - params[["cause1"]][["p"]])^exp(x_cause1 %*% params[["cause1"]][["betas"]])
    )]

    # Draw event times
    cause1_ind <- dat[["D"]] == 1
    dat[D == 1, time := generate_direct_times(runif(.N), x_cause1[cause1_ind, ], params[["cause1"]])]
    dat[D == 2, time := rweibull_KM(
      n = .N,
      shape = params[["cause2"]][["base_shape"]],
      rate = params[["cause2"]][["base_rate"]] * exp(x_cause2[!cause1_ind, ] %*% params[["cause2"]][["betas"]])
    )]

  } else if (model_type == "both_fg") {

    # Get probabilities at t -> inf
    p1 <- 1 - (1 - params[["cause1"]][["p"]])^exp(x_cause1 %*% params[["cause1"]][["betas"]])
    p2 <- 1 - (1 - params[["cause2"]][["p"]])^exp(x_cause2 %*% params[["cause2"]][["betas"]])
    if (max(p1 + p2) > 1) stop("Total failure prob exceeds once for some people eventually!")

    # Draw event indicator
    U_ind <- runif(n)
    dat[, D := fcase(
      U_ind <= p1, 1,
      p1 < U_ind & U_ind <= p1 + p2, 2,
      U_ind > p1 + p2, 0
    )]

    # Draw event times
    cause1_ind <- dat[["D"]] == 1
    cause2_ind <- dat[["D"]] == 2
    dat[D == 1, time := generate_direct_times(runif(.N), x_cause1[cause1_ind, ], params[["cause1"]])]
    dat[D == 2, time := generate_direct_times(runif(.N), x_cause2[cause2_ind, ], params[["cause2"]])]
  }

  # Add standard and normal censoring
  admin_cens <- control[["admin_cens_time"]]
  cens_rate <- control[["cens_rate"]]
  if (!is.null(admin_cens)) dat[time >= admin_cens, ':=' (D = 0, time = admin_cens)]
  if (!is.null(cens_rate)) {
    dat[, cens := rexp(.N, rate = cens_rate)]
    dat[cens < time, ':=' (D = 0, time = cens)]
  }

  return(dat)
}


process_pre_imputing <- function(dat) {

  dat[, ':=' (
    time_star = ifelse(D == 1, time, max(time) + 1e-8), # works better than setting to super large number.. @Hein?
    D_star = as.numeric(D == 1),
    miss_ind = rbinom(.N, size = 1, prob = plogis(-0.5 + Z)),
    X_compl = X,
    H_cause1 = add_marginal_cumhaz(
      timevar = time,
      statusvar = D,
      cause = 1,
      type = "cause_spec"
    ),
    H_cause2 = add_marginal_cumhaz(
      timevar = time,
      statusvar = D,
      cause = 2,
      type = "cause_spec"
    ),
    H_subdist_cause1 = add_marginal_cumhaz(
      timevar = time,
      statusvar = D,
      cause = 1,
      type = "subdist"
    )
  )]

  dat[, ':=' (
    X = factor(ifelse(miss_ind == 1, NA_character_, X)), # NA_real_ if continuous
    D = factor(D), # for imp model
    H_modif_cause1 = add_marginal_cumhaz(
      timevar = time_star,
      statusvar = D_star,
      cause = 1,
      type = "cause_spec"
    )
  )]

  setorder(dat, "time")

  return(dat)
}
