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
  kmi_parametric <<- function(time, d, x, z, dist = "weibullPH") {
    mod <- flexsurvreg(Surv(time, d == 0) ~ x + z, dist = dist,
                       inits = c(1.25, 10, -1, -1))
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
      censoring_type = "none",
      #censoring_params = list("exponential" = "1.2 * exp((as.numeric(X) - 1) + Z)"), # like 55
      params = params
    ),
    args_missingness = list(
      mech_params = list(
        "prob_missing" = 0.4,
        "mechanism_expr" = "1.5 * Z"
      )
    )
  )

  dat[, cens_time := rweibull_KM(
    n = n,
    shape = 1.25,
    rate = 10 * exp(-(as.numeric(X_obs) - 1) - Z)
  )]
  dat[, ':=' (
    D = as.numeric(time < cens_time) * as.numeric(D),
    time = pmin(cens_time, time)
  )]

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
  predmat["Lambda_01_Z", c("Lambda_01", "Z")] <- 1
  predmat["H_cens_Z", c("H_cens", "Z")] <- 1

  # Now two variants: with and without H_cens (censoring as competing event)
  predmat_compet <- predmat
  predmat_compet["X", c("Z", "D1", "Lambda_01", "H_cens", "Lambda_01_Z", "H_cens_Z")] <- 1
  predmat["X", c("Z", "D1", "Lambda_01", "Lambda_01_Z")] <- 1

  # The methods
  meths <- make.method(dat)
  meths["V"] <- paste("~I(", expression(kmi_parametric(time, D, X, Z)),")")
  meths["H_cens"] <- paste("~I(", expression(nelsaalen(V, 1 - D1, X, Z)),")")
  meths["Lambda_01"] <- paste("~I(", expression(nelsaalen(V, D1, X, Z)),")")
  #meths["H_cens"] <- paste("~I(", expression(nelsaalen(V, 1 - D1, 1, 1)),")")
  #meths["Lambda_01"] <- paste("~I(", expression(nelsaalen(V, D1, 1, 1)),")") #nelso
  meths["Lambda_01_Z"] <- "~I(Lambda_01 * Z)"
  meths["H_cens_Z"] <- "~I(H_cens * Z)"
  meths["X"] <- "logreg"

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
      "V", "Lambda_01", "Lambda_01_Z", "X"
    )
  )
  imps_compet <- mice(
    data = data.frame(dat),
    method = meths,
    m = m,
    maxit = iters,
    predictorMatrix = predmat_compet,
    visitSequence = c(
      "V", "Lambda_01", "H_cens", "Lambda_01_Z", "H_cens_Z", "X"
    )
  )

  # Bind results
  rbind(
    cbind(
      tidy(pool(with(imps, coxph(Surv(V, D1) ~ X + Z)))),
      method = "standard"
    ),
    cbind(
      tidy(pool(with(imps_compet, coxph(Surv(V, D1) ~ X + Z)))),
      method = "cens_compet"
    )
  )
}


# First with Nels, then cum bseline
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



#Other plots
rbindlist(complete(imps_compet, "all"), idcol = ".imp") |>
  #ggplot(aes(newtimes, H_subdist_cause1)) +
  ggplot(aes(H_cens, log(H_cens / Lambda_01))) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  #coord_cartesian(xlim = c(0, 2), ylim = c(0, 2)) +
  facet_wrap(~ .imp)




