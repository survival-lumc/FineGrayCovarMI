set.seed(220241030)

library(survival)
library(data.table)
library(smcfcs)
library(mice)
library(future)
library(future.apply)
library(broom)
library(ggplot2)

one_replication <- function(n = 2000) {

  # Correlated covariates
  Z <- rnorm(n = n)
  X <- rbinom(n = n, size = 1, prob = plogis(0.5 * Z))

  # MAR on Z, approx 50% missings
  missind_X <- rbinom(n = n, size = 1, prob = plogis(Z))

  # Censoring conditional on X, approx 65% will be censored
  #cens <- rexp(n = n, rate = 0.22 * exp(X))
  cens <- rexp(n = n, rate = 0.075 * exp(0.5 * X)) # more like 40% censored, less strong dep on X

  # Generate outcome n bind data
  t_tilde <- rexp(n = n, rate = 0.1 * exp(X + Z))
  time <- pmin(t_tilde, cens)
  D <- as.numeric(t_tilde < cens)
  X_miss <- ifelse(missind_X == 1, NA_real_, X)
  dat <- data.frame(time, D, X = factor(X_miss), X_obs = X, Z, missind_X)
  dat$compev_ind <- factor(as.numeric(dat$D == 0) + 1)

  # SMC-FCS settings, no rejection sampling needed since X binary
  meths_smcfcs <- make.method(dat, defaultMethod = c("norm", "logreg", "mlogit", "podds"))
  m <- 10
  iters <- 15

  # Approach 1: usual SMC-FCS
  imps_smcfcs <- smcfcs(
    originaldata = dat,
    method = meths_smcfcs,
    smtype = "coxph",
    smformula = "Surv(time, D) ~ X + Z",
    m = m,
    numit = iters
  )

  # Approach 2: SMC-FCS with censoring as competing event
  imps_smcfcs_cens <- smcfcs(
    originaldata = dat,
    method = meths_smcfcs,
    smtype = "compet",
    smformula = list(
      "Surv(time, compev_ind == 1) ~ X + Z",
      "Surv(time, compev_ind == 2) ~ X + Z"
    ),
    m = m,
    numit = iters
  )

  # NOTE: I originally did the below, which feels like it should work but it
  # goes completely wrong!

  # imps_smcfcs_cens <- smcfcs(
  #   originaldata = dat,
  #   method = meths_smcfcs,
  #   smtype = "compet",
  #   smformula = list(
  #     "Surv(time, D == 1) ~ X + Z",
  #     "Surv(time, D == 0) ~ X + Z"
  #   ),
  #   m = m,
  #   numit = iters
  # )

  # Pool
  res_smcfcs <- lapply(imps_smcfcs$impDatasets, function(impdat) {
    formula = coxph(Surv(time, D) ~ X + Z, data = impdat)
  }) |>
    pool() |>
    tidy() |>
    cbind("method" = "smcfcs", true = c(1, 1))

  res_smcfcs_cens <- lapply(imps_smcfcs_cens$impDatasets, function(impdat) {
    formula = coxph(Surv(time, D) ~ X + Z, data = impdat)
  }) |>
    pool() |>
    tidy() |>
    cbind("method" = "smcfcs_cens", true = c(1, 1))

  rbind(res_smcfcs, res_smcfcs_cens)
}

n_reps <- 100 # this will take approx 10 mins if using 3 cores

plan(multisession, workers = 3)
sims <- future_replicate(
  n = n_reps,
  expr = one_replication(n = 2000),
  simplify = FALSE
)
plan(sequential)

# For Monte-Carlo SEs
crit <- qnorm((1 - 0.95) / 2, lower.tail = FALSE)

# Plot
rbindlist(sims, idcol = "simrep") |>
  ggplot(aes(method, estimate, fill = method)) +
  #geom_violin(trim = FALSE) +
  facet_wrap(~ term, scales = "fixed") +
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
  scale_y_continuous(breaks = seq(0.8, 1.2, by = 0.05)) +
  coord_cartesian(ylim = c(0.8, 1.2)) +
  labs(y = "Pooled estimate [Monte Carlo SE of bias]")
