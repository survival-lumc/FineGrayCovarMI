set.seed(220241030)

library(survival)
library(data.table)
library(smcfcs)
library(mice)
library(future)
library(future.apply)
library(broom)
library(ggplot2)
library(flexsurv)

n <- 25000
Z <- rnorm(n = n)
X <- rbinom(n = n, size = 1, prob = plogis(0.5 * Z))

# MAR on Z, approx 50% missings
#missind_X <- rbinom(n = n, size = 1, prob = plogis(Z))
#missind_X <- rbinom(n = n, size = 1, prob = plogis(X))

# Censoring conditional on X, approx 65% will be censored
cens <- rexp(n, rate = 0.1)
#cens <- rexp(n = n, rate = 0.15 * exp(X + Z)) # more like 40% censored, less strong dep on X

# Generate outcome n bind data
t_tilde <- rexp(n = n, rate = 0.1 * exp(X + Z))
time <- pmin(t_tilde, cens)
D <- as.numeric(t_tilde < cens)
prop.table(table(D))
#X_miss <- ifelse(missind_X == 1, NA_real_, X)
dat <- data.frame(time, D, X, Z)
#dat <- data.frame(time, D, X = factor(X_miss), X_obs = X, Z, missind_X)
dat$compev_ind <- factor(as.numeric(dat$D == 0) + 1)

mod <- flexsurvreg(Surv(time, D) ~ X + Z, data = dat, dist = "exp")

simulate(
  mod,
  nsim = 1,
  newdata = dat[dat$D == 0, ],
  start = dat[dat$D == 0, "time"]
)




# Comp risks --------------------------------------------------------------


library("targets")
library("tarchetypes")
library("future")
library("future.callr")
library("here")
source(here("packages.R"))
invisible(lapply(list.files(here("R"), full.names = TRUE), source))
set.seed(89)

params <- tar_read(true_params_correct_FG_0.15)
params$cause1$betas <- c(1, 1)

one_replication <- function(n = 2000) {


  nelsaalen <<- function(timevar,
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

  kmi_parametric <<- function(time, d, x, z, dist = "exp") {
    mod <- flexsurvreg(Surv(time, d == 0) ~ x + z, dist = dist)
    potential_cens <- simulate(
      mod,
      nsim = 1,
      newdata = data.frame(x = x[d == 2], z = z[d == 2]),
      start = time[d == 2]
    )
    V_imp <- time
    V_imp[d == 2] <- potential_cens$time_1
    return(V_imp)
  }

  # Generate data - censoring depending on X and Z
  dat <- generate_dataset(
    n = n,
    args_event_times = list(
      mechanism = "correct_FG",
      censoring_type = "exponential",
      #censoring_params = list("exponential" = "0.49 * exp((as.numeric(X) - 1) + Z)"),
      censoring_params = list("exponential" = "1.2 * exp((as.numeric(X) - 1) + Z)"), # like 55
      params = params
    ),
    args_missingness = list(
      mech_params = list(
        "prob_missing" = 0.4,
        "mechanism_expr" = "1.5 * Z"
      )
    )
  )
  #prop.table(table(dat$D))
  #0       1       2
 #0.56411 0.09378 0.34211

  # Add empties
  dat[, ':=' (
    H_subdist_cause1 = NA_real_,
    H_cens = NA_real_,
    newtimes = NA_real_,
    newevent = as.numeric(D == 1)
  )]

  # Predmats
  predmat <- make.predictorMatrix(dat)
  predmat[] <- 0
  predmat["newtimes", c("time", "D", "X", "Z")] <- 1
  predmat["H_subdist_cause1", c("newtimes", "newevent")] <- 1 #marginals anyway
  predmat["H_cens", c("newtimes", "newevent")] <- 1

  # Now two variants: with and without H_cens (censoring as competing event)
  predmat_compet <- predmat
  predmat_compet["X", c("Z", "newevent", "H_subdist_cause1", "H_cens")] <- 1
  predmat["X", c("Z", "newevent", "H_subdist_cause1")] <- 1

  # The methods
  meths <- make.method(dat)
  meths["newtimes"] <- paste("~I(", expression(kmi_parametric(time, D, X, Z)),")")

  # H cens u do need to update? Try using baseline hazard instead, and with categorical
  # Z..
  meths["H_cens"] <- paste("~I(", expression(nelsaalen(newtimes, 1 - newevent)),")")

  # H subdist less important to update??
  meths["H_subdist_cause1"] <- paste("~I(", expression(nelsaalen(newtimes, newevent)),")")
  #meths["X"] <- "rf"
  meths["X"] <- "logreg" # cart


  # Let's go
  m <- 5
  iters <- 10

  imps <- mice(
    data = data.frame(dat),
    method = meths,
    m = m,
    maxit = iters,
    predictorMatrix = predmat,
    visitSequence = c("newtimes", "H_subdist_cause1", "H_cens", "X")#,
    #printFlag = FALSE
  )
  imps_compet <- mice(
    data = data.frame(dat),
    method = meths,
    m = m,
    maxit = iters,
    predictorMatrix = predmat_compet,
    visitSequence = c("newtimes", "H_subdist_cause1", "H_cens", "X")#,
    #printFlag = FALSE
  )
  # smcfcs(
  #   originaldata = data.frame(dat),
  #   smtype = "coxph",
  #   smformula = "Surv(time, newevent) ~ X + Z",
  #   m = 1,
  #   numit = 1,
  #   method = meths
  # )

  # Bind results
  rbind(
    cbind(
      tidy(pool(with(imps, coxph(Surv(newtimes, newevent) ~ X + Z)))),
      method = "standard"
    ),
    cbind(
      tidy(pool(with(imps_compet, coxph(Surv(newtimes, newevent) ~ X + Z)))),
      method = "cens_compet"
    )
  )
}

n_reps <- 50
plan(multisession, workers = 3)
sims <- future_replicate(
  n = n_reps,
  expr = one_replication(n = 2000),
  simplify = FALSE
)
plan(sequential)

crit <- qnorm((1 - 0.95) / 2, lower.tail = FALSE)

rbindlist(sims, idcol = "simrep") |>
  ggplot(aes(method, estimate, fill = method)) +
  #geom_violin(trim = FALSE) +
  facet_wrap(~ term, scales = "fixed") +
  geom_hline(aes(yintercept = 1), linetype = "dotted") +
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
  scale_y_continuous(breaks = seq(0.7, 1.3, by = 0.1)) +
  coord_cartesian(ylim = c(0.7, 1.3)) +
  labs(y = "Pooled estimate [Monte Carlo SE of bias]")





#one_replication(n = 20000)

n = 20000


plot(imps, layout = c(2, 4))

lapply(complete(imps, "all"), function(x) {
  coxph(Surv(newtimes, newevent) ~ X + Z, data = x)
}) |>
  pool()


rbindlist(complete(imps, "all"), idcol = ".imp") |>
  #ggplot(aes(newtimes, H_subdist_cause1)) +
  ggplot(aes(newtimes, H_cens)) +
  geom_step(aes(group = .imp))






# Now with exact imputation model, we need binary Z -----------------------


library("targets")
library("tarchetypes")
library("future")
library("future.callr")
library("future.apply")
library("here")
library(flexsurv)
source(here("packages.R"))
invisible(lapply(list.files(here("R"), full.names = TRUE), source))
set.seed(84899)

params <- tar_read(true_params_correct_FG_0.15)
params$cause1$betas <- c(1, 1)

# Make Z binary
generate_covariates <- function(n, X_type = "binary") {
  dat <- data.table(id = seq_len(n), Z = rbinom(n = n, size = 1L, prob = 0.5))

  dat[, X := switch(
    X_type,
    "binary" = factor(rbinom(n = n, size = 1L, prob = plogis(Z))),
    "normal" = rnorm(n = n, mean = 0.5 * Z, sd = 1)
  )]
  return(dat)
}

library(brglm2)

one_replication <- function(n = 2000) {

  # Now this return the baseline hazard rather than marginal one
  nelsaalen <<- function(timevar,
                         statusvar,
                         x, z,
                         timefix = FALSE) {
    mod <- survival::coxph(
      Surv(time, status) ~ x + z,
      control = survival::coxph.control(timefix = timefix),
      data = cbind.data.frame("time" = timevar, "status" = statusvar, "x" = x, "z" = z)
    )
    hazard <- survival::basehaz(mod)
    idx <- match(timevar, hazard[, "time"])
    return(hazard[idx, "hazard"])
  }

  # Not proper but it works
  kmi_parametric <<- function(time, d, x, z, dist = "exp") {
    #mod <- flexsurvreg(Surv(time, d == 0) ~ x + z, dist = dist,
    #                   inits = c(1.25, 10, -1, -1))
    mod <- flexsurvreg(Surv(time, d == 0) ~ x + z, dist = dist)
    potential_cens <- simulate(
      mod,
      nsim = 1,
      newdata = data.frame(x = x[d == 2], z = z[d == 2]),
      start = time[d == 2]
    )
    V_imp <- time
    V_imp[d == 2] <- potential_cens$time_1
    return(V_imp)
  }

  # Generate data - censoring depending on X and Z
  dat <- generate_dataset(
    n = n,
    args_event_times = list(
      mechanism = "correct_FG",
      censoring_type = "exponential",
      censoring_params = list(
        "exponential" = "1.2 * exp((as.numeric(X) - 1) + Z)"
        #"exponential" = "0.5 * exp((as.numeric(X) - 1) + Z + Z^2 + (as.numeric(X) - 1)^2)" # quadratic cens model
      ), # like 55
      params = params
    ),
    args_missingness = list(
      mech_params = list(
        "prob_missing" = 0.4,
        "mechanism_expr" = "1.5 * Z"
      )
    )
  )

  # Method with brglm2

  prop.table(table(dat$D))

  #prop.table(table(dat$D))
  #    0      1      2
  #0.6165 0.1130 0.2705
  #plot(prodlim(Hist(time, D == 0) ~ X + Z, data = dat))

  # Add empties
  dat[, ':=' (
    Lambda_01 = NA_real_,
    H_cens = NA_real_,
    V = NA_real_,
    D1 = as.numeric(D == 1),
    # Interaction terms
    Lambda_01_Z = NA_real_,
    H_cens_Z = NA_real_
  )]

  # Predmats
  predmat <- make.predictorMatrix(dat)
  predmat[] <- 0
  predmat["V", c("time", "D", "X", "Z")] <- 1
  predmat["Lambda_01", c("V", "D1", "X", "Z")] <- 1
  predmat["H_cens", c("V", "D1", "X", "Z")] <- 1
  #predmat["Lambda_01_Z", c("Lambda_01", "Z")] <- 1
  #predmat["H_cens_Z", c("H_cens", "Z")] <- 1

  # Now two variants: with and without H_cens (censoring as competing event)
  predmat_compet <- predmat
  #predmat_compet["X", c("Z", "D1", "Lambda_01", "H_cens", "Lambda_01_Z", "H_cens_Z")] <- 1
  #predmat["X", c("Z", "D1", "Lambda_01", "Lambda_01_Z")] <- 1
  predmat["X", c("Z", "D1", "Lambda_01", "H_cens")] <- 1

  # The methods
  meths <- make.method(dat)
  meths["V"] <- paste("~I(", expression(kmi_parametric(time, D, X, Z)),")")
  meths["H_cens"] <- paste("~I(", expression(nelsaalen(V, 1 - D1, X, Z)),")")
  meths["Lambda_01"] <- paste("~I(", expression(nelsaalen(V, D1, X, Z)),")")
  #meths["H_cens"] <- paste("~I(", expression(nelsaalen(V, 1 - D1, 1, 1)),")")
  #meths["Lambda_01"] <- paste("~I(", expression(nelsaalen(V, D1, 1, 1)),")") #nelso
  #meths["Lambda_01_Z"] <- "~I(Lambda_01 * Z)"
  #meths["H_cens_Z"] <- "~I(H_cens * Z)"
  meths["X"] <- "rf"

  #dat[, "X" := as.integer(as.numeric(X) - 1)]

  # Let's go
  m <- 10
  iters <- 5

  imps <- mice(
    data = data.frame(dat),
    method = meths,
    m = m,
    maxit = iters,
    predictorMatrix = predmat,
    visitSequence = c(
      "V", "Lambda_01", "H_cens", "X"
    ),
    remove.collinear = FALSE,
    remove.constant = FALSE
  )
  # imps_compet <- mice(
  #   data = data.frame(dat),
  #   method = meths,
  #   m = m,
  #   maxit = iters,
  #   predictorMatrix = predmat_compet,
  #   visitSequence = c(
  #     "V", "Lambda_01", "H_cens", "Lambda_01_Z", "H_cens_Z", "X"
  #   ),
  #   remove.collinear = FALSE,
  #   remove.constant = FALSE
  # )

  # Bind results
  rbind(
    cbind(
      tidy(pool(with(imps, coxph(Surv(V, D1) ~ X + Z)))),
      method = "standard"
    )#,
    # cbind(
    #   tidy(pool(with(imps_compet, coxph(Surv(V, D1) ~ X + Z)))),
    #   method = "cens_compet"
    # )
  )
}


# First with Nels, then cum bseline
n_reps <- 30
plan(multisession, workers = 3)
sims <- future_replicate(
  n = n_reps,
  expr = one_replication(n = 2000),
  simplify = FALSE
)
plan(sequential)

crit <- qnorm((1 - 0.95) / 2, lower.tail = FALSE)

rbindlist(sims, idcol = "simrep") |>
  ggplot(aes(method, estimate, fill = method)) +
  #geom_violin(trim = FALSE) +
  facet_wrap(~ term, scales = "fixed") +
  geom_hline(aes(yintercept = 1), linetype = "dotted") +
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
  scale_y_continuous(breaks = seq(0.7, 1.3, by = 0.1)) +
  coord_cartesian(ylim = c(0.7, 1.3)) +
  labs(y = "Pooled estimate [Monte Carlo SE of bias]")



#Other plots
rbindlist(complete(imps_compet, "all"), idcol = ".imp") |>
  #ggplot(aes(newtimes, H_subdist_cause1)) +
  ggplot(aes(H_cens, log(H_cens / Lambda_01))) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  #coord_cartesian(xlim = c(0, 2), ylim = c(0, 2)) +
  facet_wrap(~ .imp)





# Try some visuals --------------------------------------------------------


params <- tar_read(true_params_correct_FG_0.15)
params$cause1$betas <- c(0.75, 0.5)

dat <- generate_dataset(
  n = 200000,
  args_event_times = list(
    mechanism = "correct_FG",
    censoring_type = "exponential",
    censoring_params = list(
      #"exponential" = "1.2 * exp((as.numeric(X) - 1) + Z)"
      "exponential" = "0.49 * exp((as.numeric(X) - 1))"
      #"exponential" = "0.5 * exp((as.numeric(X) - 1) + Z + Z^2 + (as.numeric(X) - 1)^2)" # quadratic cens model
    ), # like 55
    params = params
  ),
  args_missingness = list(
    mech_params = list(
      "prob_missing" = 0,
      "mechanism_expr" = "1.5 * Z"
    )
  )
)

prop.table(table(dat$D))
dat[, V := time * (D == 1) + cens_time * (D != 1)]
dat[, D1 := factor(as.numeric(D == 1))]
dat[, F_01 := params$cause1$p * (
  1 - exp(-params$cause1$base_rate * V^params$cause1$base_shape)
)]
dat[, Lambda_01 := -log(1 - F_01)]
dat[, H_cens := 0.49 * V] # baseline

melt.data.table(
  dat,
  measure.vars = c("Lambda_01", "H_cens"),
  variable.name = "cumhaz_type",
  value.name = "cumhaz"
) |>
  ggplot(aes(cumhaz, as.numeric(X) - 1,
             col = D1, linetype = factor(Z))) +
  #geom_point(aes(shape = factor(Z))) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, bs = "cs"),
    method.args = list(family = "quasibinomial"),
    se = FALSE
  ) +
  coord_trans(y = "qlogit", xlim = c(0, 2.5)) +
  facet_wrap(~ cumhaz_type)


library(scales)

#logit_trans()
qlogit_trans <- function() trans_new(
  "qlogit", transform = function(x) qlogis(x), inverse = function(x) plogis(x)
)

dat |>
  ggplot(aes(Lambda_01, as.numeric(X) - 1,
             col = D1, linetype = factor(Z))) +
  #geom_point(aes(shape = factor(Z))) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, bs = "cs"),
    method.args = list(family = "quasibinomial"),
    se = FALSE
  ) +
  coord_trans(y = "qlogit") +
  facet_wrap(~ D1)





library(ggeffects)
mod <- glm(X ~ D1 + Z * Lambda_01, data = dat, family = binomial())
predict_response(mod, terms = c("Lambda_01", "D1", "Z")) |>
  plot()


newdf <- expand.grid(
  Z = c(0, 1),
  D1 = factor(c(0, 1)),
  Lambda_01 = seq(0, 0.20, by = 0.01)
)

cbind.data.frame(
  newdf, preds = predict(mod, newdata = newdf)
) |>
  ggplot(aes(Lambda_01, preds, col = D1)) +
  geom_point(aes(shape = factor(Z)))


tn <- trans_new



# THE sims ----------------------------------------------------------------



library("targets")
library("tarchetypes")
library("future")
library("future.callr")
library("future.apply")
library("here")
library(flexsurv)
library(brglm2)

source(here("packages.R"))
invisible(lapply(list.files(here("R"), full.names = TRUE), source))
set.seed(82249)

params <- tar_read(true_params_correct_FG_0.15)
params$cause1$betas <- c(1, 1)

# Make Z binary
generate_dataset <- function(n,
                             args_event_times,
                             args_missingness,
                             args_covariates = list(X_type = "normal")) {

  dat <- do.call(generate_covariates, args = c(list(n = n), args_covariates))
  do.call(add_event_times, args = c(list(dat = dat), args_event_times))
  do.call(add_missingness, args = c(list(dat = dat), args_missingness))
  add_cumhaz_to_dat(dat)
  return(dat[])
}

one_replication <- function(n = 2000) {

  dat <- generate_dataset(
    n = n,
    args_event_times = list(
      mechanism = "correct_FG",
      censoring_type = "exponential",
      censoring_params = list(
        "exponential" = "1.5 * exp(-X + Z + Z^2)"
        #"exponential" = "0.5 * exp((as.numeric(X) - 1) + Z + Z^2 + (as.numeric(X) - 1)^2)" # quadratic cens model
      ), # like 55
      params = params
    ),
    args_missingness = list(
      mech_params = list(
        "prob_missing" = 0.4,
        "mechanism_expr" = "1.5 * Z"
      )
    )
  )

  #prop.table(table(dat$D))
  #      0       1       2
  #0.58795 0.18906 0.22299

  # Known censoring times
  dat[, ':=' (
    V = time * (D == 1) + cens_time * (D != 1),
    D1 = as.numeric(D == 1),
    compev_ind = as.numeric(D != 1) + 1
  )]

  # SMC-FCS settings, no rejection sampling needed since X binary
  meths_smcfcs <- make.method(dat, defaultMethod = c("norm", "logreg", "mlogit", "podds"))
  m <- 5
  iters <- 20

  # Approach 1: usual SMC-FCS
  imps_smcfcs <- smcfcs(
    originaldata = dat,
    method = meths_smcfcs,
    smtype = "coxph",
    smformula = "Surv(V, D1) ~ X + Z",
    m = m,
    numit = iters,
    rjlimit = 5000
  )

  # Approach 2: SMC-FCS with censoring as competing event
  imps_smcfcs_cens <- smcfcs(
    originaldata = dat,
    method = meths_smcfcs,
    smtype = "compet",
    smformula = list(
      "Surv(V, compev_ind == 1) ~ X + Z",
      "Surv(V, compev_ind == 2) ~ X + Z + I(Z^2)"
    ),
    m = m,
    numit = iters,
    rjlimit = 5000
  )

  #plot(imps_smcfcs)
  #plot(imps_smcfcs_cens)

  # Pool
  res_smcfcs <- lapply(imps_smcfcs$impDatasets, function(impdat) {
    formula = coxph(Surv(V, D1) ~ X + Z, data = impdat)
  }) |>
    pool() |>
    tidy() |>
    cbind("method" = "smcfcs", true = c(1, 1))

  res_smcfcs_cens <- lapply(imps_smcfcs_cens$impDatasets, function(impdat) {
    formula = coxph(Surv(V, D1) ~ X + Z, data = impdat)
  }) |>
    pool() |>
    tidy() |>
    cbind("method" = "smcfcs_cens", true = c(1, 1))

  rbind(res_smcfcs, res_smcfcs_cens)
}


n_reps <- 50 # this will take approx 10 mins if using 3 cores

plan(multisession, workers = 3)
sims <- future_replicate(
  n = n_reps,
  expr = one_replication(n = 2000),
  simplify = FALSE
)
plan(sequential)


crit <- qnorm((1 - 0.95) / 2, lower.tail = FALSE)

rbindlist(sims, idcol = "simrep") |>
  ggplot(aes(method, estimate, fill = method)) +
  #geom_violin(trim = FALSE) +
  facet_wrap(~ term, scales = "free") +
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
  scale_y_continuous(breaks = seq(0.7, 1.3, by = 0.1)) +
  coord_cartesian(ylim = c(0.7, 1.3)) +
  labs(y = "Pooled estimate [Monte Carlo SE of bias]")
