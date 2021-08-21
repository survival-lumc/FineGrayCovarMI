library(riskRegression)
library(survival)
library(prodlim)
library(mstate)
library(flexsurv)
library(rstpm2)
library(ggplot2)
library(data.table)


set.seed(1202)

# Sample size, single covariate
n <- 5000#2500
X <- rbinom(n, 1, 0.5)
betaX_sd_cause1 <- log(2)
betaX_cs_cause2 <- log(1.25)

# (According to Haller)
# Define proportional subdist cause 1 and cause-spec cause 2
sd_cause1 <- Vectorize(function(t, x) {
  0.001 * exp(-(0.001 * t / log(1.5)) + betaX_sd_cause1 * x)
})

cs_cause2 <- Vectorize(function(t, x) {
  0.001 * exp(betaX_cs_cause2 * x)
})

# Function for cumulative hazard sd
cumul_sd_cause1 <- Vectorize(function(t, x) {
  log(1.5) * exp(betaX_sd_cause1 * x) * (1 - exp(-0.001 * t / log(1.5)))
})

# First define function vectorising integration in denominator
vector_denom <- Vectorize(function(t, x) {
  integrate(
    f = function(t, x) {
      sd_cause1(t, x) * exp(-cumul_sd_cause1(t, x) + cs_cause2(t, x) * t)
    },
    lower = 0,
    upper = t,
    x = x
  )$value
})

# Derive hazard of cause 2
cs_cause1 <- Vectorize(function(t, x) {

  # Compute
  num <- sd_cause1(t, x) * exp(-cumul_sd_cause1(t, x) + cs_cause2(t, x) * t)
  denom <- 1 - vector_denom(t, x)
  return(num / denom)
})

# Plot them
t_grid <- seq(1, 6000, by = 5)
plot(
  t_grid,
  sd_cause1(t_grid, x = 0),
  type = "l", col = "blue",
  ylim = c(0, 0.002),
  ylab = "Hazard",
  xlab = "Time"
)
lines(t_grid, sd_cause1(t_grid, x = 1), col = "blue", lty = "dashed")
lines(t_grid, cs_cause1(t_grid, x = 0))
lines(t_grid, cs_cause1(t_grid, x = 1), lty = "dashed")
legend(
  "topright",
  col = c("blue", "blue", "black", "black"),
  lty = rep(c("solid", "dashed"), 2),
  legend = c(
    "SD (X = 0)",
    "SD (X = 1)",
    "CS (X = 0)",
    "CS (X = 1)"
  ), bty = "n"
)

# Check cumincs
plot(
  t_grid,
  1 - exp(-cumul_sd_cause1(t_grid, x = 1)),
  type = "l",
  xlab = "Time",
  ylab = "Cumulative incidence"
)
lines(t_grid, 1 - exp(-cumul_sd_cause1(t_grid, x = 0)), lty = 2)
legend(
  "bottomright",
  col = c("black", "black"),
  lty = 1:2,
  legend = c(
    "Cause 1 (X = 1)",
    "Cause 1 (X = 0)"
  ), bty = "n"
)

# Need this to generate event times
cumul_cs_cause1 <- Vectorize(function(t, x) {
  -log(1 - vector_denom(t, x))
})


# Generate event times ----------------------------------------------------


S_min_u <- function(t, x, U) {
  exp(-cumul_cs_cause1(t, x) - cs_cause2(t, x) * t) - U
}

times_all <- rstpm2::vuniroot(
  f = S_min_u,
  U = runif(n = n),
  x = X,
  interval = cbind(rep(.Machine$double.eps, n), rep(1e5, n)),
  extendInt = "yes"
)$root

haz_cause2 <- cs_cause2(times_all, X)
delta <- 1 + rbinom(
  n,
  size = 1,
  prob = haz_cause2 / (cs_cause1(times_all, x = X) + haz_cause2)
)


# Put in data and run some modelos!
dat <- data.frame("time" = times_all, "delta" = delta, "X" = X)

# Cause-specific cause2
mod_cs_cause2 <- coxph(Surv(time, delta == 2) ~ X, data = dat)
betaX_cs_cause2
coef(mod_cs_cause2)

# FG model
mod_fg_cause1 <- FGR(Hist(time, delta) ~ X, data = dat, cause = 1)
betaX_sd_cause1
mod_fg_cause1$crrFit$coef

# Check the time-dependent effect of cause-specific 1
mod_cs_cause1 <- coxph(Surv(time, delta == 1) ~ X, data = dat)
zph_obj <- cox.zph(mod_cs_cause1, terms = TRUE)
logHR_betaX_cause1 <- log(cs_cause1(zph_obj$time, x = 1) / cs_cause1(zph_obj$time, x = 0))
coef(mod_cs_cause1)
plot(zph_obj)
lines(log(zph_obj$time), logHR_betaX_cause1, lwd = 2)
points(zph_obj$time, log(zph_obj$y))

library(magrittr)
data.frame("time" = zph_obj$time, "resid" = drop(zph_obj$y)) %>%
  ggplot(aes(time, resid)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  geom_line(aes(zph_obj$time, logHR_betaX_cause1))


# Trying non-param --------------------------------------------------------


cuminc_all <- cmprsk::cuminc(dat$time, dat$delta, group = dat$X)
plot(cuminc_all)

#cumhaz_nonparam <- -log(1 - df_unique$est)
plot(cuminc_all$`0 1`$time, -log(1 - cuminc_all$`0 1`$est), type = "l", ylim = c(0, 1), lty = "dotted")
lines(cuminc_all$`1 1`$time, -log(1 - cuminc_all$`1 1`$est), type = "l", ylim = c(0, 1), lty = "dotted")

lines(cuminc_all$`0 1`$time, cumul_sd_cause1(t = cuminc_all$`0 1`$time, x = 0))
lines(cuminc_all$`0 1`$time, cumul_sd_cause1(t = cuminc_all$`0 1`$time, x = 1))
points(dato$time, dato$sd_cumhaz_x, ylim = c(0,1))


# Try a cumul hazard
dato <- data.table(dat)
setorder(dato, time)
dato[, "risk_set" := .N - cumsum(delta == 1)]
dato[, "sd_cumhaz" := cumsum(as.numeric(delta == 1) / risk_set)]

setorder(dato, X, time)
dato[, "risk_set_x" := .N - cumsum(delta == 1), by = X]
dato[, "sd_cumhaz_x" := cumsum(as.numeric(delta == 1) / risk_set_x), , by = X]

points(dato$time, dato$sd_cumhaz_x, ylim = c(0,1))

# Try estimating back baseline hazard -------------------------------------


dat_long <- crprep(
  Tstop = "time",
  status = "delta",
  data = dat,
  trans = c(1, 2),
  keep = "X"
)

mod_long <- coxph(
  Surv(Tstart, Tstop, status == 1) ~ X,
  data = dat_long,
  subset = (failcode == 1),
  weights = weight.cens
)

# Does not work!
mod_long_param <- survreg(
  Surv(Tstart, Tstop, status == 1) ~ X,
  data = dat_long,
  subset = (failcode == 1),
  weights = weight.cens,
  dist = "exponential"
)

coef(mod_long)
mod_fg_cause1


mod_long_param <- flexsurv::flexsurvreg(
  Surv(Tstart, Tstop, status == 1) ~ X,
  data = dat_long,
  subset = (failcode == 1),
  weights = weight.cens,
  dist = "weibullPH"
)

shape <- mod_long_param$coefficients[["shape"]]
scale <- exp(-mod_long_param$coefficients[["scale"]])

exp(scale)^(-shape)


# Try rstpm2
mod_long_stpm <- rstpm2::stpm2(
  Surv(Tstart, Tstop, status == 1) ~ X,
  data = dat_long,
  subset = (failcode == 1),
  weights = weight.cens
)
