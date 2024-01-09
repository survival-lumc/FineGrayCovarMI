# Large data checks mice subdist with and without kmi?
# also with/without weighted estimator?
invisible(lapply(list.files(here("R"), full.names = TRUE), source))

failure_time_model <- "correct_FG"
args_event_times = list(
  mechanism = failure_time_model,
  censoring_type = "none", #none
  params = list(
    "cause1" = list(
      "formula" = ~ X + Z,
      "betas" = c(1, 1),
      "p" = 0.65,
      "base_rate" = 1,
      "base_shape" = 0.75
    ),
    "cause2" = list(
      "formula" = ~ X + Z,
      "betas" = c(0.75, 0.5),
      "base_rate" = 1,
      "base_shape" = 0.75
    )
  )
)

args_missingness = list(mech_params = list("prob_missing" = 0.4, "mechanism_expr" = "1.5 * Z"))
args_imputations = list(m = 10, iters = 1, rjlimit = 1000, rhs_kmi = "1")
args_covariates = list("X_type" = "binary")

set.seed(87564)
n <- 10000L
dat <- generate_dataset(
  n = n,
  args_event_times = args_event_times,
  args_missingness = args_missingness,
  args_covariates = args_covariates
)

max_ev1_time <- max(dat[D == 1]$time)
eps <- 0.1
dat[, ':=' (
  newtimes = ifelse(D == 2, max_ev1_time + eps, time),
  newevent = as.numeric(D == 1)
)]


dat[, ':=' (
  H_subdist_V = compute_marginal_cumhaz(
    timevar = newtimes,
    statusvar = newevent,
    cause = 1,
    type = "cause_spec"
  ),
  H_subdist_weighted = compute_marginal_cumhaz(
    timevar = time,
    statusvar = D,
    cause = 1,
    type = "subdist"
  )
)]

dat |>
  ggplot(aes(newtimes, H_subdist_V)) +
  geom_step(linewidth = 2, col = "red") +
  geom_step(aes(time, H_subdist_weighted), col = "blue") +
  labs(x = "Time")

dat[D == 2] |>
  ggplot(aes(H_subdist_weighted)) +
  geom_histogram(col = "black", fill = "lightblue") +
  geom_vline(xintercept = max(dat[D == 2]$H_subdist_V), col = "black", linewidth = 2) +
  coord_cartesian(xlim = c(0, ))

meths <- make.method(dat, defaultMethod = c("norm", "logreg", "polyreg", "polr"))
predmat <- make.predictorMatrix(dat)
predmat[] <- 0
predmat_V <- predmat_weighted <- predmat
predmat_V["X", c("Z", "newevent", "H_subdist_V")] <- 1
predmat_weighted["X", c("Z", "newevent", "H_subdist_weighted")] <- 1

mice_weighted <- mice(
  data = data.frame(dat),
  maxit = 1,
  m = 10,
  method = meths,
  predictorMatrix = predmat_weighted,
  printFlag = FALSE
)

mice_V <- mice(
  data = data.frame(dat),
  maxit = 1,
  m = 10,
  method = meths,
  predictorMatrix = predmat_V,
  printFlag = FALSE
)

with(mice_weighted, coxph(Surv(newtimes, newevent) ~ X + Z)) |>
  pool() |>
  tidy()

with(mice_V, coxph(Surv(newtimes, newevent) ~ X + Z)) |>
  pool() |>
  tidy()
