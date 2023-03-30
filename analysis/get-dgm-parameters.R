library(data.table)
library(riskRegression)
library(survival)
library(mice)
library(prodlim)
library(kmi)
library(smcfcs)
library(tidyverse)
library(mstate)
library(broom)
library(cmprsk)
library(flexsurv)


source("R/data-generation.R")

# Parameters
params_direct <- list(
  "cause1" = list(
    "betas" = c(0.75, 0.5),
    "p" = 0.25,
    "base_rate" = 1,
    "base_shape" = 1
  ),
  "cause2" = list(
    "betas" = c(0.75, 0.5),
    "base_rate" = 1,
    "base_shape" = 1
  )
)


large_dat <- generate_complete_dataset(
  n = 100000,
  params = params_direct,
  model_type = "direct",
  X_type = "binary",
  predictor_formulas = list("cause1" = ~ X + Z, "cause2" = ~ X + Z),
  control = list("admin_cens_time" = NULL, "cens_rate" = 1e-10) # 20% of obs are censored
)

prop.table(table(large_dat$D))
mod_fg1 <- FGR(Hist(time, D) ~ X + Z, data = large_dat, cause = 1)
mod_fg1$crrFit$coef

#mod_cox_cs1 <- coxph(Surv(time, D == 1) ~ X + Z, data = large_dat)
#mod_cox_cs1$coefficients

# Cause-specific 1 almost exactly exponential?
mod_weib_cs1 <- survreg(Surv(time, D == 1) ~ X + Z, data = large_dat)
shape_cs1 <- 1 / mod_weib_cs1$scale
base_rate_cs1 <- exp(mod_weib_cs1$coefficients[[1]])^(-shape_cs1)
lfps_cs1 <- -mod_weib_cs1$coefficients[-1] * shape_cs1

mod_weib_cs2 <- survreg(Surv(time, D == 2) ~ X + Z, data = large_dat)
shape_cs2 <- 1 / mod_weib_cs2$scale
base_rate_cs2 <- exp(mod_weib_cs2$coefficients[[1]])^(-shape_cs2)
lfps_cs2 <- -mod_weib_cs2$coefficients[-1] * shape_cs2




# Try simming new data ----------------------------------------------------

n <- 100000
Z <- rnorm(n)
X <- rbinom(n, size = 1, prob = plogis(Z))
times_cs1 <- rweibull_KM(
  n = n,
  shape = shape_cs1,
  rate = base_rate_cs1 * exp(drop(cbind(X, Z) %*% lfps_cs1))
)
times_cs2 <- rweibull_KM(
  n = n,
  shape = shape_cs2,
  rate = base_rate_cs2 * exp(drop(cbind(X, Z) %*% lfps_cs2))
)
dat_test <- cbind.data.frame(
  time = pmin(times_cs1, times_cs2),
  D = as.numeric(times_cs2 < times_cs1) + 1L,
  factor(X), Z
)

# Event proportions almost identical?
table(dat_test$D) |> prop.table()
table(large_dat$D) |> prop.table()

# Takes 33 minutes on sample size of 2022.09
system.time({
  mod_lfp <- FGR(Hist(time, D) ~ X + Z, data = dat_test, cause = 1)
  coefs_lfp <- mod_lfp$crrFit$coef
})


# beta_X ~ around 0.73, beta_Z = 0.37
print(coefs_lfp)

# No censoring!!
dt_dat <- data.table(dat_test)
dt_dat[, subdist_time := ifelse(D == 1, time, max(time) + 0.5)]
coxph(Surv(subdist_time, D == 1) ~ X + Z, data = dt_dat)


?kmi

testos_kmi <- kmi(
  Surv(time, D != 0) ~ 1,
  data = dat_test,
  etype = D,
  failcode = 1,
  nimp = m
)



# Tests with flexsurvreg --------------------------------------------------

flexsurv_cs1 <- flexsurvspline(Surv(time, D == 1) ~ X + Z, data = large_dat, k = 2)
times_cs1 <- simulate(flexsurv_cs1, nsim = 1)
flexsurv_cs2 <- flexsurvspline(Surv(time, D == 2) ~ X + Z, data = large_dat, k = 3)
times_cs2 <- simulate(flexsurv_cs2, nsim = 1)
dat_test <- cbind.data.frame(
  time = pmin(times_cs1$time_1, times_cs2$time_1),
  D = as.numeric(times_cs2$time_1 < times_cs1$time_1) + 1L,
  X = large_dat$X,
  Z = large_dat$Z
)

mod_lfp <- FGR(Hist(time, D) ~ X + Z, data = dat_test, cause = 1)
mod_lfp$crrFit$coef



test <- flexsurvreg(Surv(time, D == 1) ~ X + Z, data = large_dat, dist = "exponential")
test <- flexsurvspline(Surv(time, D == 1) ~ X + Z, data = large_dat, k = 2)

plot(
  test, type = "cumhaz",
  newdata = data.frame("X" = factor(0, levels = c("0", "1")), "Z" = 0),
  B = 0
)
plot(
  test,
  type = "hazard",
  newdata = data.frame("X" = factor(0, levels = c("0", "1")), "Z" = 0),
  B = 0
)

# Try rstpm2 --------------------------------------------------------------



library(rstpm2)
fit <- stpm2(Surv(time, D == 1) ~ X + Z, data = large_dat, df = 3)
summary(fit)

plot(
  fit,
  newdata = data.frame("X" = factor(0, levels = c("0", "1")), "Z" = 0),
  type = "haz"
)
plot(
  fit,
  newdata = data.frame("X" = factor(0, levels = c("0", "1")), "Z" = 0),
  type = "cumhaz"
)

fit_cox <- coxph(Surv(time, D == 1) ~ X + Z, data = large_dat)
cox.zph(fit_cox)
lines(basehaz(fit_cox)$time, basehaz(fit_cox)$hazard)
plot(basehaz(fit_cox))
