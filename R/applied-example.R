process_applied_dat <- function(dat_raw) {

  # Edit a few variables
  dat_raw[, ':=' (
    intdiagallo_decades = intdiagtr_allo1 / 10, # Now in decades
    year_allo1_decades = (year_allo1 - 10) / 10, # Most recent year is zero
    pb_allo1_prop = pb_allo1 / 100,
    log_hb_allo1 = log(hb_allo1) - log(10), # centered at Hb = 10
    log_wbc_allo1 = log(wbc_allo1 + 0.1) - log(15.1) # reference is WBC = 15 (use 25?)
  )]

  sm_predictors <- c(
    # These are predictors for original paper
    "hctci_risk",
    "age_allo1_decades", # Already centered at age = 60 (median age 58)
    "wbc_allo1",
    "hb_allo1",
    "pb_allo1",
    "sweat_allo1",
    "WEIGLOSS_allo1",
    "KARNOFSK_threecat",
    "PATSEX",
    "donrel_bin",
    "cmv_match",
    "ruxo_preallo1",
    "ric_allo1",
    # Add also extra donor info + auxiliary vars
    "agedonor_allo1_decades",
    "vchromos_preallo1",
    "DONSEX_allo1_1",
    "year_allo1", # year and intdiag not to be shown in forest plots
    "intdiagtr_allo1"
    #submps_allo1, tbi_allo1, source_allo1 not included
  )

  # Take subset with complete outcome data
  dat_sub <- dat_raw[, c("time_ci_adm", "status_ci_adm", sm_predictors), with = FALSE]
  dat_sub <- dat_sub[complete.cases(status_ci_adm), ]
  dat_sub[, id := seq_len(.N)]

  # Set as ordered for imputation part
  dat_sub[, ':=' (
    KARNOFSK_threecat = as.ordered(KARNOFSK_threecat),
    hctci_risk = as.ordered(hctci_risk)
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
  list("dat" = dat_sub, "sm_predictors" = sort(sm_predictors))
}


one_imputation_applied_dat <- function(dat_processed,
                                       imp_settings,
                                       ...) { # for the bootstrap tings

  options(contrasts = rep("contr.treatment", 2))

  # Load in (processed) data, make year factor for imputation of censoring times
  dat <- data.frame(dat_processed$dat)
  dat$year_allo1 <- as.factor(dat$year_allo1)
  cause <- imp_settings$cause # the outcome

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

  # Set year back to continuous in imputation model; set to decades so on same scale as others?
  dat_to_impute$year_allo1 <- as.numeric(as.character(dat_to_impute$year_allo1))

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
    cbind(dat_to_impute, "method" = "Compl. cases"),
    cbind(complete(imps_mice_subdist, "all")[[1]], "method" = "MICE subdist"),
    cbind(complete(imps_mice_csh, "all")[[1]], "method" = "MICE cause-spec"),
    cbind(imps_smcfcs_csh$impDatasets[[1]], "method" = "SMC-FCS cause-spec"),
    cbind(imps_smcfcs_subdist$impDatasets[[1]], "method" = "SMC-FCS Fine-Gray")
  )

  return(data.table(impdats))
}
