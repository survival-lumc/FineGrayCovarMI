# Extract coefficients + baseline cumulative incidence at desired time points
extract_mod_essentials <- function(fit,
                                   timepoints) {

  # For a coxph object, predictCox() gives baseline by default
  base_cuminc <- 1 - predictCox(fit, times = timepoints, centered = FALSE)$survival
  # When timepoint > max time in data, predictCox gives NA, just give max cuminc instead
  if (anyNA(base_cuminc)) base_cuminc[is.na(base_cuminc)] <- max(base_cuminc, na.rm = TRUE)
  coefs <- fit$coefficients

  essentials <- data.table(
    "time" = timepoints,
    "base_cuminc" = base_cuminc,
    "coefs" = list(coefs)
  )

  return(essentials)
}


# One replication of simulation study
one_replication <- function(args_event_times,
                            args_missingness,
                            args_imputations,
                            args_predictions,
                            args_covariates = list("X_type" = "binary"),
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
    args_missingness = args_missingness,
    args_covariates = args_covariates
  )

  # Add predictors needed for imputation
  add_cumhaz_to_dat(dat)

  # To speed up computations: assume censoring time to be known when fitting model in each imputed dataset
  # (we are comparing imputation methods, so this is okay)
  # - We do so consistently with kmi names
  if (args_event_times$censoring_type == "none") {
    max_ev1_time <- dat[D == 1, .(time = max(time))][["time"]]
    eps <- 0.1
    dat[, ':=' (
      newtimes = ifelse(D == 2, max_ev1_time + eps, time),
      newevent = as.numeric(D == 1)
    )]
  } else {
    dat[, ':=' (
      newtimes = ifelse(D == 2, cens_time, time),
      newevent = as.numeric(D == 1)
    )]
  }

  # If the censoring times are assumed to be known in the imputation phase
  # (curvy_uniform scenarios), we re-estimate subdist cumulative hazard with these
  cens_time_known <- if (args_event_times$censoring_type == "curvy_uniform") TRUE else FALSE

  if (cens_time_known) {
    dat[, H_subdist_cause1 := compute_marginal_cumhaz(
      timevar = newtimes,
      statusvar = newevent,
      cause = 1,
      type = "cause_spec" # Since we do not need weights (censoring times known)
    )]
  }

  # Analysis model formula; coxph() can be used directly
  model_formula <- update(form_rhs, Surv(newtimes, newevent) ~ .)

  # Make wrapper for function call, to not set x = TRUE each time (needed for predictCox())
  model_fun <- function(...) survival::coxph(..., x = TRUE)

  # Method 0: Full data before any missing data
  mod_full <- model_fun(update(model_formula, . ~ X_obs + . - X), data = dat)

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
  predmat_subdist["X", c("Z", "newevent", "H_subdist_cause1")] <- 1

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

  # When subdistribution time is known for the imputation,
  # .. we can use smcfcs() directly with Surv(newtimes, newevent)
  smcfcs_finegray <- if (args_event_times$censoring_type != "exponential") {
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
    # In the other case, we use kmi() + smcfcs() combination!
    df <- data.frame(dat)
    df <- df[, !(names(df) %in% c("newevent", "newtimes"))] # Avoid issues in subsetting within smcfcs()
    meths_smcfcs <- make.method(df, defaultMethod = c("norm", "logreg", "mlogit", "podds"))

    smcfcs.finegray(
      originaldata = df,
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
    method = c("mice_comp", "mice_subdist", "smcfcs_comp", "smcfcs_finegray"),
    imp_dats = c(
      list(complete(mice_comp, action = "all")),
      list(complete(mice_subdist, action = "all")),
      list(smcfcs_comp$impDatasets),
      list(smcfcs_finegray$impDatasets)
    )
  )

  # Prediction horizons read-in
  pred_times <- args_predictions$timepoints

  # Fit models in each imputed dataset - this is normally the bottleneck when model_fun was FGR()
  nested_impdats[, mods := .(
    list(lapply(imp_dats[[1]], function(imp_dat) model_fun(model_formula, data = imp_dat)))
  ), by = method]

  # Create summaries (= pooled coefficients, and the prediction 'essentials')
  summaries_impdats <- nested_impdats[, .(
    coefs_summary = list(tidy(pool(mods[[1]]), conf.int = TRUE)),
    preds_summary = list(
     rbindlist(
       lapply(mods[[1]], extract_mod_essentials, timepoints = pred_times),
       idcol = "imp"
     )
    )
  ), by = method]

  # Remove imputation objects to clear memory
  rm(mice_comp, mice_subdist, smcfcs_finegray, smcfcs_comp, nested_impdats)

  # Add summaries of other methods (CCA, and full dataset)
  method_summaries <- rbind(
    data.table(
      method = "full",
      coefs_summary = list(tidy(mod_full, conf.int = TRUE)),
      preds_summary = list(extract_mod_essentials(mod_full, pred_times))
    ),
    data.table(
      method = "CCA",
      coefs_summary = list(tidy(mod_CCA, conf.int = TRUE)),
      preds_summary = list(extract_mod_essentials(mod_CCA, pred_times))
    ),
    summaries_impdats
  )

  # Keep essential columns (remove MI-specific cols), and bind together with true betas
  # ("true" are the least-false parameters in the misspecified scenarios)
  essential_cols <- colnames(tidy(mod_CCA, conf.int = TRUE))
  method_summaries[, coefs_summary := .(
    list(cbind(coefs_summary[[1]][, essential_cols], "true" = extra_args$true_betas))
  ), by = method]

  # Eventually check if object too heavy due to nesting
  return(method_summaries)
}

# Function later to be used for pooling predicted cumulative incidences
inv_cloglog <- function(x) 1 - exp(-exp(x))
cloglog <- function(x) log(-log(1 - x))
