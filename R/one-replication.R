

# Add helper functions here --

# NOTE: we need different fun for admin censoring, which is assumed to be known!!!

# Extract coefficients + baseline cumulative incidence at desired timepoints
extract_FGR_essentials <- function(FGR_fit,
                                   timepoints) {

  # This still works even if there are factors, and you do not know their baseline level
  predictors <- all.vars(delete.response(FGR_fit$terms))
  newdat_baseline <- data.frame(matrix(data = 0L, ncol = length(predictors)))
  colnames(newdat_baseline) <- predictors
  base_cuminc <- drop(predict(FGR_fit, newdata = newdat_baseline, times = timepoints))

  # Note could ask for crr object? Not very heavy, only 10% of total FGR size
  essentials <- data.table(
    "time" = timepoints,
    "base_cuminc" = base_cuminc,
    "coefs" = list(FGR_fit$crrFit$coef)
  )

  #lapply(list(c(1, 1), c(1, 2), c(0.5, 1)), function(x) {
  #  essentials[, .(drop(unlist(coefs) %*% x)), by = c("time", "base_cuminc")]
  #})

  return(essentials)
}

# Start here
# Argument scenario ids? or pred times??
one_replication <- function(args_event_times,
                            args_missingness,
                            args_imputations,
                            args_predictions,
                            ...) {

  #
  extra_args <- list(...)

  # Generate data
  n <- 2000L
  dat <- generate_dataset(
    n = n,
    args_event_times = args_event_times,
    args_missingness = args_missingness
  )

  # If admin cens (check 5.3.3 beyersmann)
  # - change time for D = 2 to: cens_time
  # - change D to: 1 + as.numeric(D != 1)

  # Add predictors needed for imputation
  add_cumhaz_to_dat(dat)
  # Convert here once for all as df?

  model_formula <- Hist(time, D) ~ X + Z

  # Method 0: Full data
  mod_full <- FGR(Hist(time, D) ~ X_obs + Z, data = dat, cause = 1)

  # Method 0.5: Missing indicator
  #mod_missing <- FGR(Hist(time, D) ~ X + Z, data = dat, cause = 1)
  # Need to make indicator later..

  # Method 1: CCA
  mod_CCA <- FGR(model_formula, data = dat, cause = 1)

  # Method 2: MICE as if comp risks
  meths <- make.method(dat, defaultMethod = c("norm", "logreg", "polyreg", "polr"))
  predmat <- make.predictorMatrix(dat)
  predmat[] <- 0
  predmat["X", c("Z", "D", "H_cause1", "H_cause2")] <- 1

  mice_comp <- mice(
    data = data.frame(dat),
    method = meths,
    m = args_imputations$m,
    maxit = 1, #monotone
    predictorMatrix = predmat,
    printFlag = FALSE
  )

  # Method 3: SMC-FC as if comp risks
  meths_smcfcs <- make.method(dat, defaultMethod = c("norm", "logreg", "mlogit", "podds"))

  smcfcs_comp <- smcfcs(
    originaldata = data.frame(dat),
    smtype = "compet",
    smformula = list(
      "Surv(time, D == 1) ~ X + Z",
      "Surv(time, D == 2) ~ X + Z"
    ),
    method = meths_smcfcs,
    rjlimit = args_imputations$rjlimit,
    numit = args_imputations$iters,
    m = args_imputations$m
  )

  # Method 4: MICE w/ cumulative subdistribution hazard
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
  df_smcfcs_finegray <- data.frame(dat)

  browser()

  # Make indicator numeric
  df_smcfcs_finegray$D <- as.numeric(as.character(df_smcfcs_finegray$D))

  smcfcs_finegray <- smcfcs.finegray(
    originaldata = df_smcfcs_finegray,
    smformula = "Surv(time, D) ~ X + Z",
    method = meths_smcfcs, #make.method(dat, defaultMethod = c("norm", "logreg", "mlogit", "podds")),
    cause = 1,
    m = args_imputations$m,
    numit = args_imputations$iters,
    rjlimit = args_imputations$rjlimit,
  )

  ## ADD AN IFELSE HERE..
  # If cens time is known, process the data, and just just smcfcs function with method
  # coxph..

  # Create nested df with imputed datasets
  nested_impdats <- data.table(
    method = c("mice_comp", "mice_subdist", "smcfcs_finegray", "smcfcs_comp"),
    imp_dats = c(
      list(complete(mice_comp, action = "all")),
      list(complete(mice_subdist, action = "all")),
      list(smcfcs_finegray$impDatasets),
      list(smcfcs_comp$impDatasets)
    )
  )

  # Grid of months for predictions
  #pred_times <- round(seq(0, 10, by = 1 / 12), digits = 3)
  pred_times <- args_predictions$timepoints #1:10 # or every 6 months??

  # Fit models in each imputed dataset
  nested_impdats[, mods := .(
    list(lapply(imp_dats[[1]], function(imp_dat) FGR(model_formula, data = imp_dat, cause = 1)))
  ), by = method]

  # Summaries
  summaries_impdats <- nested_impdats[, .(
    coefs_summary = list(tidy(pool(lapply(mods[[1]], "[[", "crrFit")), conf.int = TRUE)),
    preds_summary = list(
      rbindlist(
        lapply(mods[[1]], extract_FGR_essentials, timepoints = pred_times),
        idcol = "imp"
      )
    ) # Imputation-specific coefs + base cuminc
  ), by = method]

  # Remove imps objects to clear memory
  rm(mice_comp, mice_subdist, smcfcs_finegray, smcfcs_comp, nested_impdats)

  # Add summaries of other methods (CCA, and full dataset)
  method_summaries <- rbind(
    data.table(
      method = "full",
      coefs_summary = list(tidy(mod_full$crrFit, conf.int = TRUE)),
      preds_summary = list(extract_FGR_essentials(mod_full, pred_times))
    ),
    data.table(
      method = "CCA",
      coefs_summary = list(tidy(mod_CCA$crrFit, conf.int = TRUE)),
      preds_summary = list(extract_FGR_essentials(mod_CCA, pred_times))
    ),
    summaries_impdats
  )

  # Keep essential columns (remove MI-specific cols), and bind together with true betas
  # ("true" are the least-false parameters in the misspecified scenarios)
  essential_cols <- colnames(tidy(mod_CCA$crrFit, conf.int = TRUE))
  method_summaries[, coefs_summary := .(
    list(cbind(coefs_summary[[1]][, essential_cols], "true" = extra_args$true_betas))
  ), by = method]

  # Add some scenario identifiers (maybs remove this? taken care by targets?)
  #method_summaries[, scen_summary := list(args_event_times)]

  # Eventually check if too heavy..
  return(method_summaries)
}


inv_cloglog <- function(x) 1 - exp(-exp(x))
cloglog <- function(x) log(-log(1 - x))
