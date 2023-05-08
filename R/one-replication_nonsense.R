# To speeeed up
one_replication_nonsense <- function(args_event_times,
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

  # Are censoring times assumed to be known
  cens_time_known <- if (
    args_event_times$censoring_type %in% c("curvy_uniform", "exponential")
  ) TRUE else FALSE

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
  #mod_full <- model_fun(update(model_formula, . ~ X_obs + Z), data = dat)

  # Method 0.5: Missing indicator
  # mod_missing <- FGR(Hist(time, D) ~ X + Z, data = dat, cause = 1)
  # (To-do at later stage)

  # Method 1: CCA , just for the tidy bit later
  mod_CCA <- model_fun(model_formula, data = dat)

  # Let's try some variants of the comp:
  # - Misspec_smcfcs_comp: same, but only X in imp model
  # - marginal smcfcs comp: no predictors?
  # - Cens as extra comp event

  dat[, "EFS_ind" := as.numeric(D != 0)]

  # Make a methods vector
  meths_smcfcs <- make.method(dat, defaultMethod = c("norm", "logreg", "mlogit", "podds"))

  # Smcfcs with even more misspec imp model
  smcfcs_misspec <- smcfcs(
    originaldata = data.frame(dat),
    smtype = "compet",
    smformula = list(
      "Surv(time, D == 1) ~ X",
      "Surv(time, D == 2) ~ X"
    ),
    method = meths_smcfcs,
    rjlimit = args_imputations$rjlimit,
    numit = args_imputations$iters,
    m = args_imputations$m
  )

  # Using the censoring as third competing event
  smcfcs_model_cens <- smcfcs(
    originaldata = data.frame(dat),
    smtype = "compet",
    smformula = list(
      "Surv(time, D == 0) ~ X + Z",
      "Surv(time, D == 1) ~ X + Z",
      "Surv(time, D == 2) ~ X + Z"
    ),
    method = meths_smcfcs,
    rjlimit = args_imputations$rjlimit,
    numit = args_imputations$iters,
    m = args_imputations$m
  )

  # SMC-FCS with the composite event (to check if it is the EFS that makes the diff)
  smcfcs_EFS <- smcfcs(
    originaldata = data.frame(dat),
    smtype = "coxph",
    smformula = "Surv(time, EFS_ind) ~ X + Z",
    method = meths_smcfcs,
    rjlimit = args_imputations$rjlimit,
    numit = args_imputations$iters,
    m = args_imputations$m
  )

  # Could even make it more nonsense with linear model of time, but lets keep
  # it respectable..

  # Create nested df with imputed datasets
  nested_impdats <- data.table(
    method = c("smcfcs_comp_omitZ", "smcfcs_cens_as_comp", "smcfcs_composite"),
    imp_dats = c(
      list(smcfcs_misspec$impDatasets),
      list(smcfcs_model_cens$impDatasets),
      list(smcfcs_EFS$impDatasets)
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

  # Remove imputation objects to clear memory
  rm(smcfcs_EFS, smcfcs_misspec, smcfcs_model_cens, nested_impdats)

  # Add summaries of other methods (CCA, and full dataset)
  method_summaries <- summaries_impdats

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
