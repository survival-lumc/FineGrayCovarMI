process_applied_dat <- function(dat_raw) {

  # Set-up predictors and auxiliary
  sm_predictors <- c(
    "age_allo1_decades", # this one is centered
    "PATSEX",
    "hctci_risk",
    "KARNOFSK_threecat",
    "ruxo_preallo1",
    "donrel_bin",
    "cmv_match",
    "wbc_allo1",
    "hb_allo1",
    "ric_allo1",
    "WEIGLOSS_allo1"
  )

  aux_predictors <- c(
    "year_allo1",
    "intdiagtr_allo1",
    "tbi_allo1"
  )

  # Subset of relevant vars
  dat_sub <- dat_raw[, c("time_ci_adm", "status_ci_adm", sm_predictors, aux_predictors), with = FALSE]

  # Other possible disease vars:
  # pb_allo1, sweat_allo1
  # Include weight loss, but avoid pb since percentage..

  # Use complete outcome data, and add id column
  dat_sub <- dat_sub[complete.cases(status_ci_adm), ]
  dat_sub[, id := seq_len(.N)]

  # Edit some variables pre-imputation
  dat_sub[, ':=' (
    KARNOFSK_threecat = as.ordered(KARNOFSK_threecat),
    hctci_risk = as.ordered(hctci_risk),
    wbc_allo1 = log(wbc_allo1 + 1) # since very skewed
  )]

  # Add the cumulative hazards
  dat_sub[, ':=' (
    H1 = compute_marginal_cumhaz(
      type = "cause_spec",
      timevar = time_ci_adm,
      statusvar = status_ci_adm,
      cause = 1
    ),
    H2 = compute_marginal_cumhaz(
      type = "cause_spec",
      timevar = time_ci_adm,
      statusvar = status_ci_adm,
      cause = 2
    ),
    # Add separate indicators
    D1 = as.numeric(status_ci_adm == 1),
    D2 = as.numeric(status_ci_adm == 2)
  )]

  # Return list of essentials
  list(
    "dat" = dat_sub,
    "sm_predictors" = sort(sm_predictors),
    "aux_predictors" = aux_predictors # think of removing this
  )
}


one_imputation_applied_dat <- function(dat_processed,
                                       imp_settings,
                                       ...) { # for the bootstrap tings

  options(contrasts = rep("contr.treatment", 2))

  # Load in (processed) data, make year factor for imputation of censoring times
  dat <- data.frame(dat_processed$dat)
  dat$year_allo1 <- as.factor(dat$year_allo1)
  cause <- imp_settings$cause # non-relapse mortality is the outcome

  # Impute the censoring times
  cens_imp <- kmi(
    formula = as.formula(paste0("Surv(time_ci_adm, status_ci_adm != 0) ~ ", imp_settings$rhs_cens)),
    data = dat,
    etype = status_ci_adm,
    failcode = cause,
    nimp = 1,
    ...
  )

  # This is the data for which covariates will be imputed
  dat_to_impute <- cbind.data.frame(cens_imp$original.data, cens_imp$imputed.data[[1]])
  dat_to_impute$newevent <- as.numeric(dat_to_impute$newevent == cause) # to avoid issues later

  # Add marginal subdistribution hazard
  dat_to_impute$H1_subdist <- compute_marginal_cumhaz(
    type = "cause_spec",
    timevar = dat_to_impute$newtimes,
    statusvar = dat_to_impute$newevent,
    cause = 1 # as a result of newevent
  )

  # Substantive model formula (for smcfcs)
  sm_predictors <- dat_processed$sm_predictors
  sm_form <- reformulate(
    response = "Surv(time_ci_adm, status_ci_adm)",
    termlabels = sm_predictors
  )

  # Set-up covariate impuation models
  meths_smcfcs <- make.method(dat_to_impute, defaultMethod = c("norm", "logreg", "mlogit", "podds"))
  meths_mice <- make.method(dat_to_impute) # pmm?

  # Predictor matrices for mice
  predmat_csh <- predmat_subdist <- make.predictorMatrix(dat_to_impute)
  predmat_csh[, setdiff(colnames(predmat_csh), c(sm_predictors, "D1", "D2", "H1", "H2"))] <- 0
  predmat_subdist[, setdiff(colnames(predmat_subdist), c(sm_predictors, "newevent", "H1_subdist"))] <- 0

  # Set year back to continuous in imputation model
  dat_to_impute$year_allo1 <- as.numeric(dat_to_impute$year_allo1) - 1

  # Run the imputations
  imps_mice_csh <- mice(
    data = dat_to_impute,
    m = 1,
    maxit = imp_settings$num_cycles,
    method = meths_mice,
    predictorMatrix = predmat_csh,
    printFlag = FALSE
  )

  imps_mice_subdist <- mice(
    data = dat_to_impute,
    m = 1,
    maxit = imp_settings$num_cycles,
    method = meths_mice,
    predictorMatrix = predmat_subdist,
    printFlag = FALSE
  )

  imps_smcfcs_csh <- smcfcs(
    originaldata = dat_to_impute,
    smtype = "compet",
    smformula = list(
      deparse1(update(sm_form, Surv(time_ci_adm, status_ci_adm == 1) ~ .)),
      deparse1(update(sm_form, Surv(time_ci_adm, status_ci_adm == 2) ~ .))
    ),
    m = 1,
    numit = imp_settings$num_cycles,
    method = meths_smcfcs,
    rjlimit = imp_settings$rjlimit
  )

  imps_smcfcs_subdist <- smcfcs(
    originaldata = dat_to_impute,
    smtype = "coxph",
    smformula = deparse1(update(sm_form, Surv(newtimes, newevent) ~ .)),
    m = 1,
    numit = imp_settings$num_cycles,
    method = meths_smcfcs,
    rjlimit = imp_settings$rjlimit
  )

  # Bind all imputed datasets, could technically nest this
  impdats <- rbind(
    cbind(complete(imps_mice_subdist, "all")[[1]], "method" = "MI subdist"),
    cbind(complete(imps_mice_csh, "all")[[1]], "method" = "MICE CSH"),
    cbind(imps_smcfcs_csh$impDatasets[[1]], "method" = "SMC-FCS CSH"),
    cbind(imps_smcfcs_subdist$impDatasets[[1]], "method" = "SMC-FCS subdist")
  )

  return(data.table(impdats))
}
