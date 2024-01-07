# Extract coefficients + baseline cumulative incidence at desired time points
extract_mod_essentials <- function(fit,
                                   timepoints) {

  # For a coxph object, predictCox() gives baseline by default
  base_cuminc <- 1 - predictCox(fit, times = timepoints, centered = FALSE)$survival

  # With this, we can predict for any reference patient and at of the timepoints
  essentials <- data.table(
    "time" = timepoints,
    "base_cuminc" = base_cuminc,
    "coefs" = list(fit$coefficients)
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

  # Whether of not censoring times are known determines much of the rest
  cens_times_known <- if (args_event_times$censoring_type == "exponential") FALSE else TRUE

  # If censoring times are known: make the V and I(D == 1) variables directly
  # (consistent with {kmi} variable naming)
  if (args_event_times$censoring_type == "none") {
    max_ev1_time <- max(dat[D == 1]$time)
    eps <- 0.1
    dat[, ':=' (
      newtimes = ifelse(D == 2, max_ev1_time + eps, time),
      newevent = as.numeric(D == 1)
    )]
  } else if (args_event_times$censoring_type == "admin") {
    dat[, ':=' (
      newtimes = ifelse(D == 2, cens_time, time),
      newevent = as.numeric(D == 1)
    )]
  }

  # If censoring times are unknown: imputed censoring times are used to estimate
  # the FG models for all methods
  obj <- if (!cens_times_known) {

    # Impute the censoring times
    kmi_obj <- do.call(
      kmi,
      args = list(
        formula = reformulate(response = "Surv(time, D != 0)", termlabels = args_imputations$rhs_kmi),
        failcode = 1,
        data = data.frame(dat),
        etype = as.symbol("D"),
        nimp = args_imputations$m
      )
    )

    # Estimate marginal cumulative subdistribution hazard in each imputed dataset
    impdats_kmi <- lapply(kmi_obj$imputed.data, function(imp_times) {
      impdat <- cbind(kmi_obj$original.data, imp_times)
      impdat$H_subdist_cause1 <- compute_marginal_cumhaz(
        timevar = impdat$newtimes,
        statusvar = impdat$newevent,
        cause = 1,
        type = "cause_spec"
      )
      impdat$newevent <- as.numeric(as.character(impdat$newevent))
      impdat
    })
    impdats_kmi
  } else {
    # When censoring times are known, can compute it directly
    dat[, H_subdist_cause1 := compute_marginal_cumhaz(
      timevar = newtimes,
      statusvar = newevent,
      cause = 1,
      type = "cause_spec"
    )]
    data.frame(dat)
  }

  # Analysis model formula: coxph() can be used directly rather than crr()
  model_formula <- update(form_rhs, Surv(newtimes, newevent) ~ .)

  # Run + summarise the (covariate) imputations
  summaries_imps <- run_imp_methods(
    obj = obj,
    model_formula = model_formula,
    args_imputations = args_imputations,
    args_predictions = args_predictions
  )

  # Run + summarise the CCA and full methods
  summaries_CCA_full <- run_CCA_full(
    obj = obj,
    model_formula = model_formula,
    args_predictions = args_predictions
  )

  # Combine the two, and add 'true' coefficient values
  method_summaries <- rbind(summaries_CCA_full, summaries_imps)
  method_summaries[, coefs_summary := .(
    list(cbind(coefs_summary[[1]], "true" = extra_args$true_betas))
  ), by = method]

  return(method_summaries)
}


# Wrapper for coxph() to not set x = TRUE each time, needed for predictCox()
model_fun <- function(...) survival::coxph(..., x = TRUE)


# Run the MI methods for the simulation study
run_imp_methods <- function(obj,
                            model_formula,
                            args_imputations,
                            args_predictions) {

  # obj is a list when the censoring times are unknown (and have been imputed with kmi)
  # .. otherwise it is a dataframe
  cens_times_known <- if (is.data.frame(obj)) TRUE else FALSE
  pred_times <- args_predictions$timepoints

  # Use the first imputed dataset to set-up imputation methods when censoring times unknown
  dat <- if (cens_times_known) obj else obj[[1]]

  # Method: MICE with cause-specific hazards
  meths <- make.method(dat, defaultMethod = c("norm", "logreg", "polyreg", "polr"))
  predmat <- make.predictorMatrix(dat)
  predmat[] <- 0
  predmat["X", c("Z", "D", "H_cause1", "H_cause2")] <- 1

  # Method: MICE with cumulative subdistribution hazard
  predmat_subdist <- predmat
  predmat_subdist[] <- 0
  predmat_subdist["X", c("Z", "newevent", "H_subdist_cause1")] <- 1

  # Method: SMC-FCS with cause-specific hazards
  meths_smcfcs <- make.method(dat, defaultMethod = c("norm", "logreg", "mlogit", "podds"))

  # Set-up general arguments for both mice() and smcfcs()
  args_smcfcs <- list(
    "method" = meths_smcfcs,
    "rjlimit" = args_imputations$rjlimit,
    "numit" = args_imputations$iters,
    "m" = args_imputations$m
  )
  args_mice <- list(
    "maxit" = 1,
    "m" = args_imputations$m,
    "printFlag" = FALSE,
    "method" = meths
  )

  # When censoring times are unknown, then we impute covariates in each kmi dataset
  if (!cens_times_known) args_mice$m <- args_smcfcs$m <- 1

  # Arguments specific to different methods
  args_mice_comp <- c(args_mice, list("predictorMatrix" = predmat))
  args_mice_subdist <- c(args_mice, list("predictorMatrix" = predmat_subdist))
  args_smcfcs_finegray <- c(
    args_smcfcs,
    list("smtype" = "coxph", "smformula" = deparse1(model_formula))
  )
  args_smcfcs_comp <- c(
    args_smcfcs,
    list(
      "smtype" = "compet",
      "smformula" = list(
        deparse1(update(model_formula, Surv(time, D == 1) ~ .)),
        deparse1(update(model_formula, Surv(time, D == 2) ~ .))
      )
    )
  )

  # Run the imputation methods, store them in nested data.table
  nested_impdats <- if (cens_times_known) {
    mice_comp <- do.call(mice, args = c(list("data" = obj), args_mice_comp))
    mice_subdist <- do.call(mice, args = c(list("data" = obj), args_mice_subdist))
    smcfcs_comp <- do.call(smcfcs, args = c(list("originaldata" = obj), args_smcfcs_comp))
    smcfcs_finegray <- do.call(smcfcs, args = c(list("originaldata" = obj), args_smcfcs_finegray))

    data.table(
      method = c("mice_comp", "mice_subdist", "smcfcs_comp", "smcfcs_finegray"),
      impdats = c(
        list(complete(mice_comp, "all")),
        list(complete(mice_subdist, "all")),
        list(smcfcs_comp$impDatasets),
        list(smcfcs_finegray$impDatasets)
      )
    )
  } else {
    imps_ls <- lapply(obj, function(impdat) {
      mice_comp <- do.call(mice, args = c(list("data" = impdat), args_mice_comp))
      mice_subdist <- do.call(mice, args = c(list("data" = impdat), args_mice_subdist))
      smcfcs_comp <- do.call(smcfcs, args = c(list("originaldata" = impdat), args_smcfcs_comp))
      smcfcs_finegray <- do.call(smcfcs, args = c(list("originaldata" = impdat), args_smcfcs_finegray))

      data.table(
        method = c("mice_comp", "mice_subdist", "smcfcs_comp", "smcfcs_finegray"),
        impdat = c(
          list(complete(mice_comp, "all")[[1]]),
          list(complete(mice_subdist, "all")[[1]]),
          list(smcfcs_comp$impDatasets[[1]]),
          list(smcfcs_finegray$impDatasets[[1]])
        )
      )
    })
    rbindlist(imps_ls)[, .(impdats = list(impdat)), by = method]
  }

  # Summarise:
  # Fit models in each imputed dataset - this is normally the bottleneck when using crr() or FGR()
  nested_impdats[, mods := .(
    list(lapply(impdats[[1]], function(impdat) model_fun(model_formula, data = impdat)))
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

  return(summaries_impdats)
}


# Run analysis on full dataset prior to missings, and CCA
run_CCA_full <- function(obj,
                         model_formula,
                         args_predictions) {

  cens_times_known <- if (is.data.frame(obj)) TRUE else FALSE
  model_formula_full <- update(model_formula, . ~ X_obs + . - X)
  pred_times <- args_predictions$timepoints

  # Same principle as with imputation methods: when censoring times are unknown,
  # use the imputed censoring times
  mod_summaries <- if (cens_times_known) {
    mod_full <- model_fun(model_formula_full, data = obj)
    mod_CCA <- model_fun(model_formula, data = obj)
    rbind(
      data.table(
        method = "full",
        coefs_summary = list(tidy(mod_full, conf.int = TRUE)),
        preds_summary = list(extract_mod_essentials(mod_full, pred_times))
      ),
      data.table(
        method = "CCA",
        coefs_summary = list(tidy(mod_CCA, conf.int = TRUE)),
        preds_summary = list(extract_mod_essentials(mod_CCA, pred_times))
      )
    )
  } else {
    nested_impdats <- rbind(
      data.table(
        method = "full",
        mods = list(lapply(obj, function(impdat) model_fun(model_formula_full, data = impdat)))
      ),
      data.table(
        method = "CCA",
        mods = list(lapply(obj, function(impdat) model_fun(model_formula, data = impdat)))
      )
    )

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
  }
  return(mod_summaries)
}


# Function later to be used for pooling predicted cumulative incidences
inv_cloglog <- function(x) 1 - exp(-exp(x))
cloglog <- function(x) log(-log(1 - x))
