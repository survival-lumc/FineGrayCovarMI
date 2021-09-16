# Covariates
#set.seed(475878)
n <- 2000
Z <- rnorm(n, sd = 1)
X <- rbinom(n, size = 1, prob = plogis(Z))
dat <- data.table(id = seq_len(n), X, Z)

# Restrict Z
dat[, "Z" := pmin(3, pmax(Z, -3))]

# Add list for parameter values
params_cuminc_cause1 <- list(
  p = 0.15,
  base_rate = 0.3,
  base_shape = 0.6,
  beta_X = 0.5,
  beta_Z = 0.25
)

params_cuminc_cause2 <- list(
  p = 0.15,
  base_rate = 0.3,
  base_shape = 0.8,
  beta_X = 0.25,
  beta_Z = 0.25
)

# Get probabilities (or draw from multinomial instead?)
dat[, ':=' (
  U_ind = runif(n), # for generating indicator
  U = runif(n), # for generating time
  p1 = 1 - (1 - params_cuminc_cause1[["p"]])^(
    exp(params_cuminc_cause1[["beta_X"]] * X + params_cuminc_cause1[["beta_Z"]] * Z)
  ),
  p2 = 1 - (1 - params_cuminc_cause2[["p"]])^(
    exp(params_cuminc_cause2[["beta_X"]] * X + params_cuminc_cause2[["beta_Z"]] * Z)
  )
)]

# Check max prob
max(dat$p1 + dat$p2)

# Assign event
dat[, D := fcase(
  U_ind <= p1, 1,
  p1 < U_ind & U_ind <= p1 + p2, 2,
  U_ind > p1 + p2, 0
)]

table(dat$D)
mean(dat$p1 + dat$p2)

# Generate time wrapper
generate_times <- function(U, x, params) {
  p <- params[["p"]]
  HR <- exp(params[["beta_X"]] * x[["X"]] + params[["beta_Z"]] * x[["Z"]])
  val <- -log(1 - (1 - (1 - U * (1 - (1 - p)^HR))^(1 / HR)) / p)
  (val / params[["base_rate"]])^(1 / params[["base_shape"]])
}


# Generate time
dat[D == 1, time := generate_times(U, .SD, params_cuminc_cause1)]
dat[D == 2, time := generate_times(U, .SD, params_cuminc_cause2)]

# Generate censoring
dat[, cens_time := rexp(.N, rate = 0.01)]
dat[D == 0, time := cens_time]
dat[cens_time <= time, ':=' (time = cens_time, D = 0)]
dat[time > 10, ':=' (time = 10, D = 0)] # admin censoring

# Check basic stats
table(dat$D)
prop.table(table(dat$D))
hist(dat$time)
hist(dat[D == 0 & time < 10]$time)

# Run models
mod1 <- FGR(Hist(time, D) ~ X + Z, cause = 1, data = dat)
mod1$crrFit$coef
unlist(params_cuminc_cause1)[c("beta_X", "beta_Z")]

mod2 <- FGR(Hist(time, D) ~ X + Z, cause = 2, data = dat)
mod2$crrFit$coef
unlist(params_cuminc_cause2)[c("beta_X", "beta_Z")]

#plot(cmprsk::cuminc(dat$time, dat$D))

# Cumulative hazard
plot(mod1$crrFit$uftime, mod1$crrFit$bfitj)



prs <- predict(mod1, newdata = data.frame("X" = 0, "Z" = 0))
plot(mod1$crrFit$uftime, as.vector(prs), type = "l", ylim = c(0, 0.4))
lines(mod1$crrFit$uftime, params_cuminc_cause1[["p"]] * (1 - exp(
  -params_cuminc_cause1[["base_rate"]] * mod1$crrFit$uftime^(params_cuminc_cause1[["base_shape"]])
)), col = "red", lwd = 2)

prs <- predict(mod2, newdata = data.frame("X" = 0, "Z" = 0))
plot(mod2$crrFit$uftime, as.vector(prs), type = "l", ylim = c(0, 0.4))
lines(mod2$crrFit$uftime, params_cuminc_cause2[["p"]] * (1 - exp(
  -params_cuminc_cause2[["base_rate"]] * mod2$crrFit$uftime^(params_cuminc_cause2[["base_shape"]])
)), col = "red", lwd = 2)



