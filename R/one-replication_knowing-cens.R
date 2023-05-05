# To speeeed up
one_replication_cens_known <- function(args_event_times,
                                       args_missingness,
                                       args_imputations,
                                       args_predictions,
                                       ...) {

  # Simple way to add true betas to summary data.frame at end
  extra_args <- list(...)

  # Keep formula for the rest
  form_rhs <- args_event_times$params$cause1$formula

  # Generate dataset
  n <- 2000L
  dat <- generate_dataset(
    n = n,
    args_event_times = args_event_times,
    args_missingness = args_missingness
  )

  # Add predictors needed for imputation
  add_cumhaz_to_dat(dat)

  # Are censoring times assumed to be known?
  cens_time_known <- if (args_event_times$censoring_type == "curvy_uniform") TRUE else FALSE

  # If censoring times are known: set competing event times to their censoring time
  # The marginal cumulative subdistribution hazard is re-estimated based on this
  # new time variables and D == 1 indicator
  if (cens_time_known) {
    # Consistency with kmi names
    dat[, ':=' (
      newtimes = ifelse(D == 2, cens_time, time),
      newevent = D_star
    )]
    dat[, H_subdist_cause1 := compute_marginal_cumhaz(
      timevar = newtimes,
      statusvar = newevent,
      cause = 1,
      type = "cause_spec"
    )]
  }

  # Analysis model formula; coxph() is used directly when cens_time_known == TRUE
  model_formula <- if (cens_time_known) {
    update(form_rhs, Surv(newtimes, newevent) ~ .)
  } else update(form_rhs, Hist(time, D) ~ .)

  # Make wrapper for function call
  model_fun <- if (cens_time_known) {
    function(...) survival::coxph(..., x = TRUE)
  }  else {
    function(...) riskRegression::FGR(..., cause = 1)
  }

  # Method 0: Full data before any missing data
  mod_full <- model_fun(update(model_formula, . ~ X_obs + Z), data = dat)

  # Method 0.5: Missing indicator
  # mod_missing <- FGR(Hist(time, D) ~ X + Z, data = dat, cause = 1)
  # (To-do at later stage)

  # Method 1: CCA
  mod_CCA <- model_fun(model_formula, data = dat)

  # Method 2: MICE with cause-specific hazards
  meths <- make.method(dat, defaultMethod = c("norm", "logreg", "polyreg", "polr"))
  predmat <- make.predictorMatrix(dat)
  predmat[] <- 0
  predmat["X", c("Z", "D", "H_cause1", "H_cause2")] <- 1

  mice_comp <- mice(
    data = data.frame(dat),
    method = meths,
    m = args_imputations$m,
    maxit = 1, # monotone
    predictorMatrix = predmat,
    printFlag = FALSE
  )

  # Method 3: SMC-FCS with cause-specific hazards
  meths_smcfcs <- make.method(dat, defaultMethod = c("norm", "logreg", "mlogit", "podds"))

  smcfcs_comp <- smcfcs(
    originaldata = data.frame(dat),
    smtype = "compet",
    smformula = list(
      deparse1(update(form_rhs, Surv(time, D == 1) ~ .)),
      deparse1(update(form_rhs, Surv(time, D == 2) ~ .))
    ),
    method = meths_smcfcs,
    rjlimit = args_imputations$rjlimit,
    numit = args_imputations$iters,
    m = args_imputations$m
  )

  # Method 4: MICE with cumulative subdistribution hazard
  predmat_subdist <- predmat
  predmat_subdist[] <- 0
  predmat_subdist["X", c("Z", "D_star", "H_subdist_cause1")] <- 1

  mice_subdist <- mice(
    data = data.frame(dat),
    method = meths,
    m = args_imputations$m,
    maxit = 1,
    predictorMatrix = predmat_subdist,
    printFlag = FALSE
  )

  # Method 5: SMC-FCS Fine-Gray
  dat[, D := as.numeric(as.character(D))] # Make indicator numeric

  # When censoring times are known, we can use smcfcs() directly
  # .. with Surv(time_star, D_star) outcome
  smcfcs_finegray <- if (cens_time_known) {
    smcfcs(
      originaldata = data.frame(dat),
      smformula = deparse1(update(form_rhs, Surv(newtimes, newevent) ~ .)),
      method = meths_smcfcs,
      smtype = "coxph",
      m = args_imputations$m,
      numit = args_imputations$iters,
      rjlimit = args_imputations$rjlimit
    )
  } else {
    # In the other case, we use kmi() + smcfcs() combination
    smcfcs.finegray(
      originaldata = data.frame(dat),
      smformula = deparse1(update(form_rhs, Surv(time, D) ~ .)),
      method = meths_smcfcs,
      cause = 1,
      m = args_imputations$m,
      numit = args_imputations$iters,
      rjlimit = args_imputations$rjlimit
    )
  }

  # Create nested df with imputed datasets
  nested_impdats <- data.table(
    method = c("mice_comp", "mice_subdist", "smcfcs_comp"),
    imp_dats = c(
      list(complete(mice_comp, action = "all")),
      list(complete(mice_subdist, action = "all")),
      list(smcfcs_comp$impDatasets)
    )
  )

  # Prediction horizons
  pred_times <- args_predictions$timepoints

  # Fit models in each imputed dataset; this is time-consuming when model_fun = FGR()
  nested_impdats[, mods := .(
    list(lapply(imp_dats[[1]], function(imp_dat) model_fun(model_formula, data = imp_dat)))
  ), by = method]

  # Create summaries (= pooled coefficients, and the prediction 'essentials')
  summaries_impdats <- nested_impdats[, .(
    coefs_summary = list(tidy(pool_tweaked(mods[[1]]), conf.int = TRUE)),
    preds_summary = list(
      rbindlist(
        lapply(mods[[1]], extract_mod_essentials, timepoints = pred_times),
        idcol = "imp"
      )
    )
  ), by = method]

  # We need to do this separately for smcfcs finegray, since we should make
  # .. use of the imputed censoring times
  mods_smcfcs_finegray <- lapply(
    smcfcs_finegray$impDatasets,
    function(imp) coxph(update(form_rhs, Surv(newtimes, newevent) ~ .), data = imp, x = TRUE)
  )
  summaries_smcfcs_finegray <- data.table(
    "method" = "smcfcs_finegray",
    "coefs_summary" = list(tidy(pool_tweaked(mods_smcfcs_finegray), conf.int = TRUE)),
    "preds_summary" = list(
      rbindlist(
        lapply(mods_smcfcs_finegray, extract_mod_essentials, timepoints = pred_times),
        idcol = "imp"
      )
    )
  )

  # Remove imputation objects to clear memory
  rm(mice_comp, mice_subdist, smcfcs_finegray, smcfcs_comp, nested_impdats)

  # Add summaries of other methods (CCA, and full dataset)
  method_summaries <- rbind(
    data.table(
      method = "full",
      coefs_summary = list(tidy_tweaked(mod_full, conf.int = TRUE)),
      preds_summary = list(extract_mod_essentials(mod_full, pred_times))
    ),
    data.table(
      method = "CCA",
      coefs_summary = list(tidy_tweaked(mod_CCA, conf.int = TRUE)),
      preds_summary = list(extract_mod_essentials(mod_CCA, pred_times))
    ),
    summaries_impdats,
    summaries_smcfcs_finegray
  )

  # Keep essential columns (remove MI-specific cols), and bind together with true betas
  # ("true" are the least-false parameters in the misspecified scenarios)
  essential_cols <- colnames(tidy_tweaked(mod_CCA, conf.int = TRUE))
  method_summaries[, coefs_summary := .(
    list(cbind(coefs_summary[[1]][, essential_cols], "true" = extra_args$true_betas))
  ), by = method]

  # Eventually check if object too heavy due to nesting
  return(method_summaries)
}

# Function later to be used for pooling predicted cumulative incidences
inv_cloglog <- function(x) 1 - exp(-exp(x))
cloglog <- function(x) log(-log(1 - x))
