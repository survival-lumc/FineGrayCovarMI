# Libraries
library(data.table)
library(prodlim)
library(riskRegression)
library(survival)
library(mice)
library(ggplot2)


# Data-generating functions -----------------------------------------------


# Covariates
set.seed(84862)
n <- 500
Z <- rnorm(n, sd = 1)
X <- rbinom(n, size = 1, prob = plogis(Z))
dat <- data.table(id = seq_len(n), X, Z)

# Restrict
dat[, "Z" := pmin(3, pmax(Z, -3))]

# Set parameters
params_sd_cause1 <- list(
  "base_rate" = 0.15,
  "base_cuminc" = 0.25, # Baseline cumulative incidence for t -> infinity
  "beta_X" = 0.75,
  "beta_Z" = 0.25
)

params_cs_cause2 <- list(
  "base_rate" = 0.25, #0.01, #0.05,
  "base_shape" = 0.5, #0.75, # Decreasing hazard
  "gamma_X" = 0.5,
  "gamma_Z" = -0.25
)

# Make helper functions
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
  } else {
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
      # Check covar values
      covars <- paste(round(c(x[["X"]], x[["Z"]]), 2), collapse = ", ")
      stop(paste0("Likely integral itself becomes infinite. Covars: ", covars))
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


# Checks
t_grid <- 1:2600
combo <- list("X" = 1, "Z" = -3)

plot(
  t_grid,
  haz_cs_cause1(
    t_grid,
    params_sd_cause1 = params_sd_cause1,
    params_cs_cause2 = params_cs_cause2,
    x = combo
  ),
  type = "l"
)

plot(
  t_grid,
  integral_haz_cause1(
    t_grid,
    params_sd_cause1 = params_sd_cause1,
    params_cs_cause2 = params_cs_cause2,
    x = combo
  ),
  type = "l"
)

haz_cs_cause1(
  10000000,
  params_sd_cause1 = params_sd_cause1,
  params_cs_cause2 = params_cs_cause2,
  x = combo
)


# Plots -------------------------------------------------------------------


# Check the hazards
hazards_df <- data.table(
  expand.grid(
    "time_grid" = seq(1, 10, by = 1),
    "X" = c(0, 1),
    "Z" = c(-3, -2, -1, 0, 1, 2, 3)
  )
)

hazards_df[, ':=' (
  "haz_sd_cause1" = haz_sd_cause1(time_grid, list("X" = X, "Z" = Z), params_sd_cause1),
  "haz_cs_cause2" = haz_cs_cause2(time_grid, list("X" = X, "Z" = Z), params_cs_cause2),
  "haz_cs_cause1" = haz_cs_cause1(time_grid, list("X" = X, "Z" = Z), params_sd_cause1, params_cs_cause2)
), by = c("X", "Z")]

melt.data.table(
  hazards_df,
  id.vars = c("time_grid", "X", "Z"),
  value.name = "hazard",
  variable.name = "hazard_type"
) |>
  ggplot(aes(time_grid, hazard)) +
  geom_line(
    aes(col = hazard_type, group = hazard_type, linetype = hazard_type),
    size = 1.25,
    alpha = 0.8
  ) +
  coord_cartesian(
    #ylim = c(0, 0.1)
  ) +
  facet_grid(X ~ Z) +
  theme_bw()


# Generating times --------------------------------------------------------


S_min_u <- function(t, x, params_sd_cause1, params_cs_cause2, U) {
  H1 <- haz_cs_cause1(t, x, params_sd_cause1, params_cs_cause2, type = "cumulative")
  H2 <- haz_cs_cause2(t, x, params_cs_cause2, type = "cumulative")
  exp(-H1 - H2) - U
}


# Checking root of this
root_denom <- function(t, x, params_sd_cause1, params_cs_cause2) {
  1 - integral_haz_cause1(t, x, params_sd_cause1, params_cs_cause2)
}

root_denom(t = 5, list("X" = 0, "Z" = 0), params_sd_cause1, params_cs_cause2)

uniroot(
  f = root_denom,
  interval = c(.Machine$double.eps, 1e4),
  extendInt = "yes",
  params_sd_cause1 = params_sd_cause1,
  params_cs_cause2 = params_cs_cause2,
  x = list("X" = 1, "Z" = 2.5)#, check.conv =
)

generate_time <- function(params_sd_cause1,
                          params_cs_cause2,
                          x,
                          interval_upper = 1e4) {
  # Generate u
  U <- runif(1)

  # Evaluate survival at limit
  surv_upper <- S_min_u(
    interval_upper,
    params_sd_cause1 = params_sd_cause1,
    params_cs_cause2 = params_cs_cause2,
    x = x,
    U = U
  )

  if (is.nan(surv_upper)) {
    stop("This aint workin'bud.")
  } else if (surv_upper > 0) {
    interval_upper # and then censor by maxt
  } else {
    # Generate time by inversion
    uniroot(
      f = S_min_u,
      interval = c(.Machine$double.eps, interval_upper),
      extendInt = "yes",
      params_sd_cause1 = params_sd_cause1,
      params_cs_cause2 = params_cs_cause2,
      x = x,
      U = U
    )$root
  }
}


replicate(5, expr = {
  n <- 500
  Z <- rnorm(n, sd = 1)
  X <- rbinom(n, size = 1, prob = plogis(Z))
  dat <- data.table(id = seq_len(n), X, Z)

  # Restrict
  dat[, "Z" := pmin(3, pmax(Z, -3))]

  dat[, times := generate_time(
    params_sd_cause1 = params_sd_cause1,
    params_cs_cause2 = params_cs_cause2,
    x = .SD,
    interval_upper = 1000
  ), by = id]

}, simplify = F)

dat[, times := generate_time(
  params_sd_cause1 = params_sd_cause1,
  params_cs_cause2 = params_cs_cause2,
  x = .SD,
  interval_upper = 15
), by = id]

hazards_dat <- dat[, .(
  haz1 = haz_cs_cause1(times, .SD, params_sd_cause1, params_cs_cause2),
  haz2 = haz_cs_cause2(times, .SD, params_cs_cause2)
), by = id]

delta <- 1 + rbinom(
  n,
  size = 1,
  prob = hazards_dat$haz2 / (hazards_dat$haz1 + hazards_dat$haz2)
)

maxt <- 10
dat[, D := delta]
dat[times >= maxt, ':=' (
  times = maxt,
  D = 0
)]
table(dat$D)

FGR(Hist(times, D) ~ X + Z, cause = 1, data = dat)
coxph(Surv(times, D == 2) ~ X + Z, data = dat)
plot(cmprsk::cuminc(dat$times, dat$D), xlim = c(0, 10))
ggplot(dat, aes(times)) +
  geom_histogram(fill = "blue", col = "black") +
  facet_wrap(~ D, scales = "free")

dat_bis <- copy(dat)
dat_bis[, integ := integral_haz_cause1(
  times, x = .SD, params_sd_cause1, params_cs_cause2), by = id]

hist(dat_bis$integ)

cumhaz_cs_cause1 <- function(t, x, params_sd_cause1, params_cs_cause2) {
  -log(1 - integral_haz_cause1(t, x, params_sd_cause1, params_cs_cause2))
}

tapply(cbind.data.frame("X" = x, "Z" = z),
       x = c(0, 1),
       z = seq(-5, 5, by = 0.25),
       SIMPLIFY = FALSE)

# Test which covariate values are messing this up
# Error handling
# See https://stackoverflow.com/questions/21084624/data-table-and-error-handling-using-try-statement

# Monitor all functions over time; or just the integral
