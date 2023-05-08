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

  # Now the nonsense methods: lets make numeric indicators for both competing
  # events
  dat[, ':=' (
    D_ev1 = as.numeric(D == 1),
    D_ev2 = as.numeric(D == 2)
  )]
  #debug(smcfcs:::smcfcs.core)

  # Make a methods vector
  meths_smcfcs <- make.method(dat, defaultMethod = c("norm", "logreg", "mlogit", "podds"))

  # Impute as if Cox just cause 1 (ignoring competing risks)
  smcfcs_surv_ev1 <- smcfcs(
    originaldata = data.frame(dat),
    smtype = "coxph",
    smformula = "Surv(time, D_ev1) ~ X + Z",
    method = meths_smcfcs,
    rjlimit = args_imputations$rjlimit,
    numit = args_imputations$iters,
    m = args_imputations$m
  )

  # Impute as if Cox just cause 2 (!!) (ignoring competing risks, and wrong event)
  smcfcs_surv_ev2 <- smcfcs(
    originaldata = data.frame(dat),
    smtype = "coxph",
    smformula = "Surv(time, D_ev2) ~ X + Z",
    method = meths_smcfcs,
    rjlimit = args_imputations$rjlimit,
    numit = args_imputations$iters,
    m = args_imputations$m
  )

  # Could even make it more nonsense with linear model of time, but lets keep
  # it respectable..

  # Create nested df with imputed datasets
  nested_impdats <- data.table(
    method = c("smcfcs_cox_ev1", "smcfcs_cox_ev2"),
    imp_dats = c(
      list(smcfcs_surv_ev1$impDatasets),
      list(smcfcs_surv_ev2$impDatasets)
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
  rm(smcfcs_surv_ev1, smcfcs_surv_ev2, nested_impdats)

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
