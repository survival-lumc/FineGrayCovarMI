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


# Helpers for direct methods ----------------------------------------------





# Helpers for CS/indirect methods -----------------------------------------


haz_sd_cause1 <- function(t, x, params, type = "hazard") {
  base <- log(1 / (1 - params[["base_cuminc"]]))
  lp <- params[["beta_X"]] * x[["X"]] + params[["beta_Z"]] * x[["Z"]]

  if (type == "hazard") {
    params[["base_rate"]] * exp(-params[["base_rate"]] * t / base + lp)
  } else { # Cumulative
    exp(lp) * base * (1 - exp(-params[["base_rate"]] * t / base))
  }
}

haz_cs_cause2 <- function(t, x, params, type = "hazard") {
  lp <- params[["gamma_X"]] * x[["X"]] + params[["gamma_Z"]] * x[["Z"]]
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
