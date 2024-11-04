library(data.table)
library(survival)
library(mice)
library(smcfcs)
library(future)
library(future.apply)
library(broom)
library(ggplot2)

#https://onlinelibrary.wiley.com/doi/10.1002/sim.4016
set.seed(433)

rweibull_KM <- function(n, alph, lam) {
  (-log(stats::runif(n)) / lam)^(1 / alph)
}

one_replication <- function(n = 250) {

  # Independent covariates
  Zm <- rnorm(n = n, mean = 1 / 2, sd = 1 / 12)
  Zo <- rbinom(n = n, size = 1, prob = 0.5) # they are indep! Maybe more bias here
  #.. but probs because of lack of compatibility + continuous Zm

  betas <- c(-log(2), log(2))
  shape_baseline <- 0.8
  rate_baseline <- 0.6^0.8

  t_tilde <- rweibull_KM(
    n = n,
    alph = shape_baseline,
    lam = rate_baseline * exp(betas[1] * Zm + betas[2] * Zo)
  )
  mean_cens <- (0.5 * Zm + 0.5 * Zo) / 2.1 # big big effects
  cens <- rexp(n = n, rate = 1 / mean_cens)
  times <- pmin(t_tilde, cens) # = X in the paper
  delta <- as.numeric(t_tilde < cens)
  missind_Zm <- rbinom(
    n = n,
    size = 1,
    prob = plogis(0.35 - 2 * delta - times - Zo) # ~ 35% missings
  )
  Zm_amp <- ifelse(missind_Zm == 1, NA_real_, Zm)
  dat <- data.frame(times, delta, Zm = Zm_amp, Zo, Zm_orig = Zm)

  # barely 20% events, think this is riddled with mc error

  meths_smcfcs <- make.method(dat, defaultMethod = c("norm", "logreg", "mlogit", "podds"))
  imps_smc <- smcfcs(
    originaldata = dat,
    method = meths_smcfcs,
    smtype = "coxph",
    smformula = "Surv(times, delta) ~ Zm + Zo", # include missing indicator in outcome model?
    m = 10,
    numit = 15,
    rjlimit = 2500
  )
  plot(imps_smc)

  res_imps <- lapply(imps_smc$impDatasets, function(impdat) {
    coxph(Surv(times, delta) ~ Zm + Zo, data = impdat)
  }) |>
    pool() |>
    tidy() |>
    cbind(meth = "smcfcs")

  res_CCA <- tidy(coxph(Surv(times, delta) ~ Zm + Zo, data = dat)) |>
    cbind(meth = "CCA")

  cbind(rbindlist(list(res_imps, res_CCA), fill = TRUE), "true" = rep(betas, times = 2))
}

#one_replication(n = 250)



# Started 16:08
n_sim <- 1000 # 1000
plan(multisession, workers = 3)
sims <- future_replicate(
  n = n_sim,
  expr = one_replication(n = 1000),
  #expr = one_replication(n = 5000),
  simplify = FALSE
)
plan(sequential)

crit <- qnorm((1 - 0.95) / 2, lower.tail = FALSE)

rbindlist(sims, idcol = "simrep")[, rel_bias := 100 * (estimate - true) / true] |>
  ggplot(aes(meth, estimate - true, fill = meth)) +
  geom_violin(trim = FALSE) +
  facet_wrap(~ term, scales = "fixed") +
  geom_hline(aes(yintercept = 0), linetype = "dotted") +
  stat_summary(
    fun = function(x) mean(x),
    fun.min = function(x) mean(x) - crit * (sd(x) / sqrt(length(x))),
    fun.max = function(x) mean(x) + crit * (sd(x) / sqrt(length(x))),
    geom = "crossbar",
    alpha = 0.5,
    fill = "black",
    fatten = 1,
    linewidth = 0.25
  ) +
  theme_bw()# +
 # coord_cartesian(ylim = c(-60, 60)) +
  #scale_y_continuous(breaks = c(-50, -25, -10, -5, 5, 10, 25, 50))

rbindlist(sims, idcol = "simrep")[, rel_bias := 100 * (estimate - true) / true][, .(
  sd(estimate  - true) / sqrt(.N),
  sd(estimate  - true)
), by = c("meth", "term")]

ests <- rbindlist(sims, idcol = "simrep")
#ests[, cum_bias := cumsum(estimate - true) / simrep, by = c("meth", "term")]
test <- ests[,  .(
  mcse = unlist(lapply(seq_len(.N), function(i) {
    sd(estimate[seq_len(i)] - true[seq_len(i)]) / sqrt(i)
  })),
  sim_rep = seq_len(.N),
  estimate,
  true
), by = c("meth", "term")]
test[, cum_bias := cumsum(estimate - true) / sim_rep, by = c("meth", "term")]

test |>
  ggplot(aes(sim_rep, cum_bias, col = term, linetype = term)) +
  geom_line() +
  facet_wrap(~ meth) +
  theme_minimal()

test[!is.na(mcse)] |>
  ggplot(aes(sim_rep, mcse, col = term, linetype = term)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ meth) +
  theme_minimal() +
  geom_hline(yintercept = 0.01)

ests[, ':=' (
  test = mean(estimate[seq_len(simrep)])
), by = c("meth", "term")]

ests[,  .(
  lapply(seq_len(.N), function(i) cbind(mean(estimate[seq_len(i)] - true[seq_len(i)]), i))
), by = c("meth", "term")]

ests[, .SD[seq_len()], by = c("meth", "term")]

ests[, .(
  lapply(seq_len(.N), function(i) {
    list("x" = 1, "y" 1)
  })
)]


