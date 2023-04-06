source(here("R/data-generation.R"))
source(here("R/smcfcs.finegray.R"))

args_event_times <- list(
  mechanism = "correct_FG",
  params = tar_read(true_params_correct_FG),
  censoring_type = "none"
)
args_missingness <- list(mech_params = list("prob_missing" = 0.1, "mechanism_expr" = "Z"))
args_imputations <- list(m = 5, iters = 1, rjlimit = 1000) #note rjlimit does not matter for binary X


# Add helper functions here --

# Extract coefficients + baseline cumulative incidence at desired timepoints
extract_FGR_essentials <- function(FGR_fit,
                                   timepoints) {
  predictors <- all.vars(delete.response(FGR_fit$terms))
  coefs <- FGR_fit$crrFit$coef

  # This still works even if there are factors, and you do not know their baseline level
  newdat_baseline <- data.frame(matrix(data = 0L, ncol = length(predictors)))
  colnames(newdat_baseline) <- predictors
  base_cuminc <- setNames(
    drop(predict(FGR_fit, newdata = newdat_baseline, times = timepoints)),
    nm = as.character(timepoints)
  )

  # Note could ask for crr object? Not very heavy, only 10% of total FGR size
  essentials <- list(
    "coef" = coefs,
    "base_cuminc" = base_cuminc
  )
  return(essentials)
}

# Input is a nested data.table with imputed datasets
summarise_imputations <- function(nested_dat,
                                  timepoints) {

  # Fit models in each imputed dataset
  nested_dat[, mods := .(
    list(lapply(imp_dats[[1]], function(imp_dat) FGR(model_formula, data = imp_dat, cause = 1)))
  ), by = method]

  # Get essentials
  nested_dat[, .(
    list(
      lapply(mods[[1]], extract_FGR_essentials, timepoints = pred_times),
      tidy(pool(lapply(mods[[1]], "[[", "crrFit")), conf.int = TRUE)
    )
  ), by = method][, .(list(V1)), by = method]

}

# Start here
one_replication <- function(args_event_times,
                            args_missingness,
                            args_imputations) {

  # Generate data
  n <- 2000L
  dat <- generate_dataset(
    n = n,
    args_event_times = args_event_times,
    args_missingness = args_missingness
  )

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
  df_smcfcs_finegray$D <-  as.numeric(
    factor(df_smcfcs_finegray$D, levels = c(0, levels(df_smcfcs_finegray$D)))
  ) - 1L

  #df_smcfcs_finegray <- as.
  smcfcs_finegray <- smcfcs.finegray(
    originaldata = df_smcfcs_finegray,
    smformula = "Surv(time, D) ~ X + Z",
    method = meths_smcfcs, #make.method(dat, defaultMethod = c("norm", "logreg", "mlogit", "podds")),
    cause = 1,
    m = args_imputations$m,
    numit = args_imputations$iters,
    rjlimit = args_imputations$rjlimit,
  )

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

  # Remove imps objects to clear memory
  rm(mice_comp, mice_subdist, smcfcs_finegray, smcfcs_comp)

  # Grid of months for predictions
  pred_times <- round(seq(0, 10, by = 1 / 12), digits = 3)

  # Fit models in each imputed dataset
  nested_impdats[, mods := .(
    list(lapply(imp_dats[[1]], function(imp_dat) FGR(model_formula, data = imp_dat, cause = 1)))
  ), by = method]

  # Get essentials
  nested_impdats[, .(
    list(
      lapply(mods[[1]], extract_FGR_essentials, timepoints = pred_times),
      tidy(pool(lapply(mods[[1]], "[[", "crrFit")), conf.int = TRUE)
    )
  ), by = method][, .(list(V1)), by = method]



  #extract_FGR_essentials(mod_CCA, timepoints = pred_times)



  nested_impdats

  # Summaries for CCA and full
  dt_CCA <- data.table(
    "method" = "CCA",
    "mods" = list(mod_CCA),
    "summaries" = list(
      list(
        pred_components = ,
        pooled_summary = tidy(mod_CCA$crrFit, conf.int = TRUE)
      )
    )
  )

  # Or should we generate all sim dats first (for one scenario), then nest? See mcdermott blog post..
  summs <- rbindlist(
    list(
      mi_summaries,
      data.table(
        method = "CCA",
        mods = list(mod_CCA),
        coefs_summ = list(tidy(mod_CCA$crrFit, conf.int = TRUE))
      ),
      data.table(
        method = "full",
        mods = list(mod_full),
        coefs_summ = list(tidy(mod_full$crrFit, conf.int = TRUE))
      )
    ),
    fill = TRUE
  )


    # Use .SD above??

    # morisot pooling
    # https://github.com/survival-lumc/CauseSpecCovarMI/blob/master/R/illustrative-analysis-helpers.R

    summs
    summs2 <- copy(summs)
    summs2[, mods := NULL] #FGR is what takes up memory!!
    #https://cran.r-project.org/web/packages/tidyr/vignettes/nest.html
    # object size
    #https://stackoverflow.com/questions/70878796/how-to-store-nested-data-efficiently-in-r
    # https://osf.io/f6pxw/download

    # Note predicted cuminc only at event 1!! Possible efficient alternative:
    # store baseline hazard at given timepoints?? but then have to do in each imp dataset..

    # Bind the true ones..too
    # intersect() colnames for binding


  return(nested_imps)
}


test_imps <- replicate(
  n = 10,
  expr = {
    one_replication(
      args_event_times,
      args_missingness,
      args_imputations
    )
  },
  simplify = FALSE
)

df <- rbindlist(test_imps)
pryr::object_size(df)
pryr::mem_used()
# https://win-vector.com/2014/05/30/trimming-the-fat-from-glm-models-in-r/


# Use to get baseline cumulative incidence at grid of monthly timepoints
1 - (1 - predict(mod_full, newdata = list(X_obs = 0, Z = 0),
                 times = c(3, 4, 5)))^exp(mod_full$crrFit$coef[1])
