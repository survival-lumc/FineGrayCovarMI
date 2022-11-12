library(riskRegression)
library(prodlim)
library(survival)

# Covars
set.seed(9877890)
n <- 2000
X <- rbinom(n, size = 1, prob = 0.5)

# Exponential baseline hazards for all (notation Saadati)
alpha_20 <- 0.5
alpha_10 <- 0.1
h_10 <- 0.05 # Not possible to have exponential subdist hazard??
beta <- 0.5
HR <- exp(beta * X)

# Generate latent times
T1 <- rexp(n, rate = alpha_10 * HR)

# Invert cause 2
U <- runif(n)
T2 <- log(U) / (h_10 * (1 - HR) - alpha_10 * (1 - HR) - alpha_20)
time <- pmin(T1, T2)
eps <- as.numeric(T2 < T1) + 1

#
dat <- cbind.data.frame(X, time, eps)
table(dat$eps)

# Look at coefs
cs1 <- coxph(Surv(time, eps == 1) ~ X, data = dat)
fg1 <- FGR(Hist(time, eps) ~ X, data = dat, cause = 1)
cs2 <- coxph(Surv(time, eps == 2) ~ X, data = dat)

cs1$coefficients
fg1$crrFit$coef

plot(cox.zph(cs2))
cox.zph(cs2)



# Try with weibull --------------------------------------------------------


# Covars
set.seed(9890)
n <- 5000
X <- rbinom(n, size = 1, prob = 0.5)

# Exponential baseline hazards for all (notation Saadati)
alpha_20 <- 0.1
alpha_10 <- 0.1
beta <- 0.5

S_min_u <- function(t, U, x) {
  cumhaz <- t * alpha_20 + (1 - exp(beta * x)) *
    (t * alpha_10 - 0.05 * t^0.75)
  exp(-cumhaz) - U
}

T2 <- mapply(function(u, x) {
  stats::uniroot(
    S_min_u,
    interval = c(.Machine$double.eps, 1000),
    extendInt = "yes",
    U = u,
    x = x
  )$`root`
}, u = runif(n), x = X)

# Generate latent times
T1 <- rexp(n, rate = alpha_10 * exp(beta * X))

time <- pmin(T1, T2)
eps <- as.numeric(T2 < T1) + 1

#
dat <- cbind.data.frame(X, time, eps)
table(dat$eps)

# Look at coefs
cs1 <- coxph(Surv(time, eps == 1) ~ X, data = dat)
fg1 <- FGR(Hist(time, eps) ~ X, data = dat, cause = 1)
cs2 <- coxph(Surv(time, eps == 2) ~ X, data = dat)
plot(cox.zph(cs2))
cox.zph(cs2)

cs1$coefficients
fg1$crrFit$coef
