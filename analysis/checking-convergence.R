# Check convergence of smcfcs FG

source(here::here("packages.R"))
invisible(lapply(list.files(here::here("R"), full.names = TRUE), source))
options(contrasts = rep("contr.treatment", 2))

dat_processed <- tar_read(applied_dat)
dat <- data.frame(dat_processed$dat)
sm_predictors <- dat_processed$sm_predictors
sm_predictors

imp_settings <- list(
  num_imputations = 10, #100,
  num_cycles = 10, #30,
  num_batches = 10,#10,
  rjlimit = 10000,
  rhs_cens = c("year_allo1_decades + intdiagallo_decades + tceldepl_bin"),
  cause = 1 # relapse
)

coxph(
  Surv(time_ci_adm, status_ci_adm == 0) ~ year_allo1_decades +
    age_allo1_decades +
    PATSEX +
    tceldepl_bin +
    cmv_match +
    donrel_bin +
    ric_allo1 +
    submps_allo1 +
    intdiagallo_decades,
  data = dat
) |> summary()

cens_imp <- kmi(
  formula = as.formula(paste0("Surv(time_ci_adm, status_ci_adm != 0) ~ ", imp_settings$rhs_cens)),
  data = dat,
  etype = status_ci_adm,
  failcode = imp_settings$cause,
  nimp = imp_settings$num_imputations
)

library(future)
library(future.apply)
set.seed(875955)

plan("multisession", workers = 3)
imps_FG <- future_lapply(cens_imp$imputed.data, FUN = function(cens_times) {

  dat_to_impute <- cbind.data.frame(cens_imp$original.data, cens_times)
  dat_to_impute$newevent <- as.numeric(dat_to_impute$newevent == imp_settings$cause)
  meths_smcfcs <- make.method(dat_to_impute, defaultMethod = c("norm", "logreg", "mlogit", "podds"))
  dat_to_impute$year_allo1_decades <- as.numeric(as.character(dat_to_impute$year_allo1_decades))
  sm_form <- reformulate(
    response = "Surv(time_ci_adm, status_ci_adm)",
    termlabels = sm_predictors
  )

  smcfcs(
    originaldata = dat_to_impute,
    smtype = "coxph",
    smformula = deparse1(update(sm_form, Surv(newtimes, newevent) ~ .)),
    m = 1,
    numit = imp_settings$num_cycles,
    method = meths_smcfcs,
    rjlimit = imp_settings$rjlimit
  )
}, future.seed = TRUE)
plan(sequential)

imps_comb <- smcfcs:::combine_smcfcs_objects(imps_FG)
imps_comb$impDatasets |> length()
#plot(imps_comb)

lapply(imps_comb$impDatasets, function(impdat) {
  coxph(reformulate(
    response = "Surv(newtimes, newevent)",
    #response = "Surv(time_ci_adm, status_ci_adm == 2)",
    termlabels = sm_predictors
  ), data = impdat)
}) |>
  pool() |>
  tidy() |>
  ggforestplot::forestplot(
    name = term, estimate = estimate,
    se = std.error,
    logodds = TRUE,
    ds = TRUE,
    xtickbreaks = c(0.7, 0.9, 0.8, 1, 1.25, 1.5, 2),
    xlim = c(-0.5, 2.25)
  )

