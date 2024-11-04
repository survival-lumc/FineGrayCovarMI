library(data.table)
library(survival)
library(mice)
library(smcfcs)
library(future)
library(future.apply)
library(broom)
library(ggplot2)

set.seed(4148645)

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

one_replication <- function(setting = c("cens_cond_X",
                                        "cens_missind_X",
                                        "cens_missind_mixture"),
                            n = 2000) {

  # Correlated covariates
  Z <- rnorm(n = n)
  X <- rbinom(n = n, size = 1, prob = plogis(0.5 * Z))
  missind_X <- rbinom(n = n, size = 1, prob = plogis(1 + Z)) #  ~MAR-Z, 50% missings

  # Different kinds of cens; approx 70% cens
  cens <- if (setting == "cens_cond_X") {
    #rexp(n = n, rate = 0.22 * exp(1 * X))
    rexp(n = n, rate = 0.1 * exp(1 * X))
  } else if (setting == "cens_missind_X") {
    rexp(n = n, rate = 0.22 * exp(1 * missind_X))
  } else if (setting == "cens_missind_mixture") {
    # only depends on X for those missing X; can make more complicated mixture
    rexp(n = n, rate = 0.27 * exp(1 * missind_X * X))
  }

  # Generate outcome n bind data
  t_tilde <- rexp(n = n, rate = 0.1 * exp(0.5 * X + 0.5 * Z))
  time <- pmin(t_tilde, cens)
  D <- as.numeric(t_tilde < cens)
  X_miss <- ifelse(missind_X == 1, NA_real_, X)
  dat <- data.frame(time, D, X = factor(X_miss), X_obs = X, Z, missind_X)

  # Add some things
  dat$cens_ind <- as.numeric(dat$D == 0)
  dat$cumhaz <- nelsaalen(dat$time, dat$D)
  dat$cumhaz_cens <- nelsaalen(dat$time, dat$cens_ind) #cens ind is 1 - D, so too correlated

  # Just the usual imp with smcfcs
  #meths_smcfcs <- make.method(dat, defaultMethod = c("norm", "logreg", "mlogit", "podds"))
  meths <- make.method(dat, defaultMethod = c("rf", "logreg", "polyreg", "polr")) # or rf?
  predmat <- predmat_cens <- make.predictorMatrix(dat)
  predmat[] <- predmat_cens[] <- 0
  predmat["X", c("Z", "D", "cumhaz")] <- 1
  predmat_cens["X", c("Z", "D", "cumhaz", "cumhaz_cens")] <- 1

  imps_mice <- mice(
    data = data.frame(dat),
    method = meths,
    m = 10,
    maxit = 1,
    predictorMatrix = predmat,
    printFlag = FALSE
  )
  imps_mice_cens <- mice(
    data = data.frame(dat),
    method = meths,
    m = 10,
    maxit = 1,
    predictorMatrix = predmat_cens,
    printFlag = FALSE
  )

  #approach <- "compet"
  # imps_smc <- if (approach == "standard") {
  #   imps_smc <- smcfcs(
  #     originaldata = dat,
  #     method = meths_smcfcs,
  #     smtype = "coxph",
  #     smformula = "Surv(time, D) ~ X + Z", # include missing indicator in outcome model?
  #     m = 10,
  #     numit = 15
  #   )
  # } else {
  #   smcfcs(
  #     originaldata = dat,
  #     method = meths_smcfcs,
  #     smtype = "compet",
  #     smformula = list(
  #       "Surv(time, D == 1) ~ X + Z",
  #       "Surv(time, D == 0) ~ X + Z"
  #     ),
  #     m = 10,
  #     numit = 15
  #   )
  # }

  # Bring results together
  res_mice <- with(imps_mice, coxph(Surv(time, D) ~ X + Z)) |>
    pool() |>
    tidy() |>
    cbind("method" = "imps_mice")

  res_mice_cens <- with(imps_mice_cens, coxph(Surv(time, D) ~ X + Z)) |>
    pool() |>
    tidy() |>
    cbind("method" = "imps_mice_cens")

  res_cca <- cbind(tidy(coxph(Surv(time, D) ~ X + Z, data = dat)), "method" = "CCA")
  cbind(
    rbindlist(list(res_mice, res_mice_cens, res_cca), fill = TRUE),
    "true" = rep(c(.5, .5), times = 3)
  )
}


plan(multisession, workers = 3)
sims <- future_replicate(
 n = 200,
 expr = one_replication(setting = c("cens_cond_X")),
 simplify = FALSE
)
plan(sequential)


crit <- qnorm((1 - 0.95) / 2, lower.tail = FALSE)

rbindlist(sims) |>
  ggplot(aes(method, estimate, fill = method)) +
  geom_violin(trim = FALSE) +
  facet_wrap(term ~ .) +
  geom_hline(aes(yintercept = true), linetype = "dotted") +
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
  theme_bw()


# Mixture has bias, but tinyyyyy
# Cond ~ X some bias too, only for X
# miss ind X nothing
# .. include missind X as auxiliary var??


scens <- c("cens_cond_X",
           "cens_missind_X",
           "cens_missind_mixture")

plan(multisession, workers = 3)
sims <- future_lapply(scens, function(scen) {
  replicate(
    n = 200,
    expr = one_replication(setting = scen),
    simplify = FALSE
  ) |> rbindlist(idcol = "sim_rep") |> cbind("setting" = scen)
}, future.seed = TRUE)
plan(sequential)


#sims <- future_replicate(
#  n = 250,
#  expr = one_replication(setting = c("cens_cond_X")),
#  simplify = FALSE
#)

crit <- qnorm((1 - 0.95) / 2, lower.tail = FALSE)

rbindlist(sims) |>
  ggplot(aes(method, estimate, fill = method)) +
  geom_violin(trim = FALSE) +
  facet_grid(term ~ setting) +
  geom_hline(aes(yintercept = true), linetype = "dotted") +
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
  theme_bw() #+
  #scale_y_continuous(breaks = seq(0.75, 1.25, by = 0.05))


# Things to try:
# - with mice (because of speed: w/ or w/o the nelson aalen of cens in imp model)
# - Cens ~ X
# - Cens ~ only on X for those with miss ind
# - Cens ~ mixture of diff rates for with and without missing X (both depending ~ X)
# - Cens ~ X for those with missing data, Cens ~ 1 for those without
