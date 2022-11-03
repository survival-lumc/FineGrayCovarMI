# Libraries
library(data.table)
library(riskRegression)
library(survival)
library(mice)
library(prodlim)
library(kmi)
library(smcfcs)


# Data-generating functions -----------------------------------------------


# Covariates
set.seed(5681)
n <- 2500
Z <- rnorm(n, sd = 0.25)
X <- rbinom(n, size = 1, prob = plogis(Z))
dat <- data.table(X, Z)

# Depending on effects + range of X.. this is no reliable

# Betas
betas_sd1 <- c(0.5, 0.1)
betas_cs2 <- c(-0.5, -0.1)

# Calculate linear predictors
lp_sd1 <- drop(as.matrix(dat) %*% betas_sd1)
lp_cs2 <- drop(as.matrix(dat) %*% betas_cs2)

# Functions
sd_cause1 <- Vectorize(function(t, lp_sd1) {
  0.001 * exp(-(0.001 * t / log(1.5)) + lp_sd1)
})

cs_cause2 <- Vectorize(function(t, lp_cs2) {
  0.001 * exp(lp_cs2)
})

cumul_sd_cause1 <- Vectorize(function(t, lp_sd1) {
  log(1.5) * exp(lp_sd1) * (1 - exp(-0.001 * t / log(1.5)))
})

# First define function vectorising integration in denominator
vector_denom <- Vectorize(function(t, lp_sd1, lp_cs2) {
  integrate(
    f = function(t, lp_sd1, lp_cs2) {
      sd_cause1(t, lp_sd1) * exp(-cumul_sd_cause1(t, lp_sd1) + cs_cause2(t, lp_cs2) * t)
    },
    lower = 0,
    upper = t,
    lp_sd1 = lp_sd1,
    lp_cs2 = lp_cs2
  )$value
})

cs_cause1 <- Vectorize(function(t, lp_sd1, lp_cs2) {

  # Compute
  num <- sd_cause1(t, lp_sd1) * exp(-cumul_sd_cause1(t, lp_sd1) + cs_cause2(t, lp_cs2) * t)
  denom <- 1 - vector_denom(t, lp_sd1, lp_cs2)
  return(num / denom)
})

cumul_cs_cause1 <- Vectorize(function(t, lp_sd1, lp_cs2) {
  -log(1 - vector_denom(t, lp_sd1, lp_cs2))
})

S_min_u <- Vectorize(function(t, lp_sd1, lp_cs2, U) {
  exp(-cumul_cs_cause1(t, lp_sd1, lp_cs2) - cs_cause2(t, lp_cs2) * t) - U
})

times_all <- rstpm2::vuniroot(
  f = S_min_u,
  U = runif(n = n),
  lp_sd1 = lp_sd1,
  lp_cs2 = lp_cs2,
  interval = cbind(rep(.Machine$double.eps, n), rep(1e5, n)),
  extendInt = "yes"
)$root

haz_cause2 <- cs_cause2(times_all, lp_cs2)
delta <- 1 + rbinom(
  n,
  size = 1,
  prob = haz_cause2 / (cs_cause1(times_all, lp_sd1, lp_cs2) + haz_cause2)
)


combo <- c(1, 20) # c(X, Z)
ts <- seq(0.1, 10000, by = 25)
plot(
  ts,
  1 - exp(-cumul_sd_cause1(ts, drop(betas_sd1 %*% combo))),
  type = "l",
  ylab = "Cumulative incidence"
)
abline(h = 1, lty = "dotted")


# Check convergence at t = 0 of sd1 and cs1
plot(
  ts,
  sd_cause1(ts, drop(betas_sd1 %*% combo)), type = "l",
  ylab = "Hazrd"
)
abline(v = 0, lty = "dashed")
lines(
  ts,
  cs_cause2(ts,
            drop(betas_cs2 %*% combo)), lty = "dotdash"
)
lines(
  ts,
  cs_cause1(ts, drop(betas_sd1 %*% combo),
            drop(betas_cs2 %*% combo)), lty = "dotted"
)


# Generate and check parameters -------------------------------------------

dat_full <- data.table(X, Z, delta, "time" = times_all)

FGR(Hist(time, delta) ~ X + Z, data = dat_full, cause = 1)
coxph(Surv(time, delta == 2) ~ X + Z, data = dat_full)


# Induce missing data -----------------------------------------------------

dat_full[, miss_ind := rbinom(n, 1, prob = plogis(-Z))]
dat_full[, X_obs := ifelse(miss_ind == 1, NA, X)]


# Testing imputations -----------------------------------------------------

FGR(Hist(time, delta) ~ X + Z, data = dat_full, cause = 1)

nelsaalen_timefixed <- function(dat,
                                timevar,
                                statusvar,
                                timefix = FALSE) {

  timevar <- as.character(substitute(timevar))
  statusvar <- as.character(substitute(statusvar))
  time <- dat[, timevar]
  status <- dat[, statusvar]
  mod <- survival::coxph(
    Surv(time, status) ~ 1,
    control = survival::coxph.control(timefix = timefix)
  )
  hazard <- survival::basehaz(mod)
  idx <- match(time, hazard[, "time"])
  return(hazard[idx, "hazard"])
}

# Make indicators
dat_full[, ':=' (
  ind_ev1 = as.numeric(delta == 1),
  ind_ev2 = as.numeric(delta == 2),
  ind_not1 = as.numeric(delta != 1)
)]

dat_full[, ':=' (
  H1 = nelsaalen_timefixed(data.frame(.SD), "time", "ind_ev1"),
  H2 = nelsaalen_timefixed(data.frame(.SD), "time", "ind_ev2")
)]

setorder(dat_full, time)
dat_full[, "risk_set" := .N - cumsum(delta == 1)]
dat_full[, "SD1" := cumsum(as.numeric(delta == 1) / risk_set)]
dat_full[, "risk_set" := NULL]

# Make interactions
dat_full[, ':=' (
  SD1_ev1 = ind_ev1 * SD1,
  H2_not1 = ind_not1 * H2
)]

# Initial mice
predmat <- matrix(
  0, ncol(dat_full), ncol(dat_full),
  dimnames = list(colnames(dat_full), colnames(dat_full))
)
predmat["X", c("Z")] <- 1
