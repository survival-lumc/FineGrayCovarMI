set.seed(2241030)

library(survival)
library(data.table)
library(smcfcs)
library(mice)
library(future)
library(future.apply)
library(broom)
library(ggplot2)
library(splines)

nelsaalen <- function(timevar,
                      statusvar,
                      timefix = FALSE) {
  mod <- survival::coxph(
    Surv(time, status) ~ 1,
    control = survival::coxph.control(timefix = timefix),
    data = cbind.data.frame("time" = timevar, "status" = statusvar)
  )
  hazard <- survival::basehaz(mod)
  idx <- match(timevar, hazard[, "time"])
  return(hazard[idx, "hazard"])
}

rweibull_KM <- function(n, shape, rate) (-log(runif(n)) / rate)^(1 / shape)

one_replication <- function(n = 2000, meth_X = "norm") {#, mtry = function(p) p) {

  # Correlated covariates
  Z <- rnorm(n = n)
  X <- rnorm(n = n, mean = 0.5 * Z) # this is probs important

  # MAR on Z, approx 50% missings
  missind_X <- rbinom(n = n, size = 1, prob = plogis(0.5 * Z))
  prop.table(table(missind_X))

  # Censoring conditional on X AND Z, approx 70% will be censored
  cens <- rweibull_KM(n = n, shape = 4, rate = 0.1 * exp(-X - Z))
  #cens <- runif(n, min = 20, max = 40)

  # Generate outcome n bind data
  #t_tilde <- rexp(n = n, rate = 0.1 * exp(X + Z))
  t_tilde <- rweibull_KM(n = n, shape = 0.75, rate = 0.12 * exp(X + Z))
  time <- pmin(t_tilde, cens)
  D <- as.numeric(t_tilde < cens)
  prop.table(table(D))
  X_miss <- ifelse(missind_X == 1, NA_real_, X)
  #dat <- data.frame(time, D, X = factor(X_miss), X_obs = X, Z, missind_X)
  dat <- data.frame(time, D, X = X_miss, X_obs = X, Z, missind_X)
  dat$compev_ind <- factor(as.numeric(dat$D == 0) + 1)

  # SMC-FCS settings, no rejection sampling needed since X binary
  #meths_smcfcs <- make.method(dat, defaultMethod = c("norm", "logreg", "mlogit", "podds"))
 # m <- 5
  #iters <- 15

  # Approach 1: usual SMC-FCS
  # imps_smcfcs <- smcfcs(
  #   originaldata = dat,
  #   method = meths_smcfcs,
  #   smtype = "coxph",
  #   smformula = "Surv(time, D) ~ X + Z",
  #   m = m,
  #   numit = iters
  # )

  # Approach 2: SMC-FCS with censoring as competing event
  # imps_smcfcs <- smcfcs(
  #   originaldata = dat,
  #   method = meths_smcfcs,
  #   smtype = "compet",
  #   smformula = list(
  #     "Surv(time, compev_ind == 1) ~ X + Z",
  #     "Surv(time, compev_ind == 2) ~ X + Z"
  #   ),
  #   m = m,
  #   numit = iters
  # )

  # Both are exponential and so proportional? Maybs choose different dists
  dat$H_ev <- nelsaalen(dat$time, dat$D)
  dat$H_cens <- nelsaalen(dat$time, 1 - dat$D)

  #formz <- lapply(make.formulas(dat), function(form) update(form, . ~ 1))
  #formz_compet <- formz
  #formz$X <- X ~ D + H_ev * Z
  #formz_compet$X <- X ~ D + H_ev * Z + H_cens * Z

  #form1 <- list(X ~ D * H_ev * Z)
  #form1 <- name.formulas(form1)
  #form2 <- list(X ~ D * H_ev * Z * H_cens)
  #form2 <- name.formulas(form2)

  predmat <- make.predictorMatrix(dat)
  predmat[] <- 0
  predmat_compet <- predmat
  predmat_compet["X", c("Z", "D", "H_ev", "H_cens")] <- 1
  predmat["X", c("Z", "D", "H_ev")] <- 1

  # The methods
  meths <- make.method(dat)
  #meths <- list("X" = "")
  meths["X"] <- meth_X # cart instead?

  m <- 10
  iters <- 1

  imps <- mice(
    data = data.frame(dat),
    method = meths,
    m = m,
    maxit = iters,
    #formulas = form1,
    predictorMatrix = predmat,
    printFlag = FALSE
  )
  imps_compet <- mice(
    data = data.frame(dat),
    method = meths,
    m = m,
    maxit = iters,
    #formulas = form2,
    predictorMatrix = predmat_compet,
    printFlag = FALSE
  )

  # Bind results
  rbind(
    cbind(
      tidy(pool(with(imps, coxph(Surv(time, D) ~ X + Z)))),
      method = "standard",
      model = "main",
      true = 1
    ),
    cbind(
      tidy(pool(with(imps_compet, coxph(Surv(time, D) ~ X + Z)))),
      method = "cens_compet",
      model = "main",
      true = 1
    ),
    cbind(
      tidy(pool(with(imps, coxph(Surv(time, D == 0) ~ X + Z)))),
      method = "standard",
      model = "cens",
      true = -1
    ),
    cbind(
      tidy(pool(with(imps_compet, coxph(Surv(time, D == 0) ~ X + Z)))),
      method = "cens_compet",
      model = "cens",
      true = -1
    )
  )
}

n_reps <- 150 # this will take approx 10 mins if using 3 cores

plan(multisession, workers = 3)
sims <- future_replicate(
  n = n_reps,
  expr = one_replication(n = 2000, meth_X = "rf"),
  simplify = FALSE
)
plan(sequential)

# For Monte-Carlo SEs
crit <- qnorm((1 - 0.95) / 2, lower.tail = FALSE)

# Plot
rbindlist(sims, idcol = "simrep") |>
  ggplot(aes(method, estimate, fill = method)) +
  #geom_violin(trim = FALSE) +
  facet_grid(model ~ term, scales = "free_y") +
  geom_hline(aes(yintercept = true), linetype = "dotted") +
  # For the monte carlo boxes
  stat_summary(
    fun = function(x) mean(x),
    fun.min = function(x) mean(x) - crit * (sd(x) / sqrt(length(x))),
    fun.max = function(x) mean(x) + crit * (sd(x) / sqrt(length(x))),
    geom = "crossbar",
    alpha = 0.5,
    aes(fill = method),
    fatten = 1,
    linewidth = 0.25
  ) +
  theme_bw() +
  #scale_y_continuous(breaks = seq(0.8, 1.2, by = 0.05)) +
  #coord_cartesian(ylim = c(0.4, 1.25)) +
  labs(y = "Pooled estimate [Monte Carlo SE of bias]")
