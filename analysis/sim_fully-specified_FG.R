n <- 5000
p <- 0.4
dat <- data.table(X = rbinom(n, size = 1, 0.5), U = runif(n))
dat[, D := 1 + rbinom(n, size = 1, prob = (1 - p)^exp(0.5 * X))]
dat[D == 2, time := rexp(.N, rate = exp(0.25 * X))]
dat[D == 1, time := -log(
  1 - (1 - (1 - U * (1 - (1 - p)^exp(0.5 * X)))^(1 / exp(0.5 * X))) / p
)]
mod <- FGR(Hist(time, D) ~ X, cause = 1, data = dat)
mod$crrFit$coef

coxph(Surv(time, D == 2) ~ X, data = dat, subset = D == 2)


# Try with two covariates -------------------------------------------------


# Helpers
rweibull_KM <- function(n, alph, lam) {
  (-log(stats::runif(n)) / lam)^(1 / alph)
}

# Sample size and limit of F1(t) at infinity
n <- 5000
p <- 0.5

# Generate covariates
dat <- data.table(Z = rnorm(n, sd = 0.25), U = runif(n))
dat[, X := rbinom(n, size = 1, prob = plogis(Z))]

dat[, D := 1 + rbinom(n, size = 1, prob = (1 - p)^exp(0.5 * X - 0.25 * Z))]
dat[, prob_d_two := (1 - p)^exp(0.5 * X - 0.25 * Z)]
dat[, prob_d_one := 1 - prob_d_two]

# Try with weibull event 2
#dat[D == 2, time := rexp(.N, rate = exp(0.25 * X - 0.1 * Z))]
dat[D == 2, time := rweibull_KM(.N, alph = 0.75, lam = exp(-0.5 * X + 0.1 * Z))]

dat[D == 1, time := -log(
  1 - (1 - (1 - U * (1 - (1 - p)^exp(0.5 * X - 0.25 * Z)))^(1 / exp(0.5 * X - 0.25 * Z))) / p
)]
mod <- FGR(Hist(time, D) ~ X + Z, cause = 1, data = dat)
mod$crrFit$coef

mod <- coxph(Surv(time, D == 2) ~ X + Z, data = dat, subset = D == 2)
coef(mod)
plot(cox.zph(mod))

# Just for fun
coxph(Surv(time, D == 1) ~ X + Z, data = dat)
coxph(Surv(time, D == 2) ~ X + Z, data = dat)

hist(dat[D == 2][["time"]])
hist(dat[D == 1][["time"]])

# Is this it???

plot(cmprsk::cuminc(dat$time, dat$D))
cumo <- cmprsk::cuminc(dat$time, dat$D)
cumo$`1 1`$est |> max()
cumo$`1 2`$est |> max()
dat$prob_d_one |> mean()
dat$prob_d_two |> mean()
dat$prob_d_one |> hist()
dat$prob_d_two |> hist()
