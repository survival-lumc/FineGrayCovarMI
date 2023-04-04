args_event_times <- list(
  mechanism = "correct_FG",
  params = tar_read(true_params_correct_FG),
  censoring_type = "none"
)
args_missingness <- list(mech_params = list("prob_missing" = 0.1, "mechanism_expr" = "Z"))
args_imputations <- list(m = 2, iters = 2, rjlimit = 1000) #note rjlimit does not matter for binary X

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
  nested_imps <- data.table(
    method = c(
      "mice_comp",
      "mice_subdist",
      "smcfcs_finegray",
      "smcfcs_comp"
    ),
    imps = c(
      list(complete(mice_comp, action = "all")),
      list(complete(mice_subdist, action = "all")),
      list(smcfcs_finegray$impDatasets),
      list(smcfcs_comp$impDatasets)
    )
  )

  # Or

  # Fit substantive model in all datasets - write wrappers for the lapplys?
  nested_imps_mods <- nested_imps[, .(
    mods = list(
      lapply(imps[[1]], function(imp_dat) FGR(model_formula, data = imp_dat, cause = 1))
    )
  ), by = method]

  # Summarise MI info: substantive models + pooled summary
  mi_summaries <- nested_imps_mods[, .(
    mods = mods,
    coefs_summ = list(
      cbind(method, tidy(pool(lapply(mods[[1]], function(FGR_fit) FGR_fit$crrFit)), conf.int = TRUE))
    )
  ), by = method]

  rbindlist(mi_summaries$coefs_summ)


  # Test if this is quicker?
  imp_dats <- list(
    "mice_comp" = complete(mice_comp, action = "all"),
    "mice_subdist" = complete(mice_subdist, action = "all"),
    "smcfcs_finegray" = smcfcs_finegray$impDatasets,
    "smcfcs_comp" = smcfcs_comp$impDatasets
  )

  rbindlist(imp_dats)

  # Do the same for CCA methods
 test <-  nested_imps[, .(imp_dats = imps[[1]]), by = method]
 test[, imp_ind := seq_len(.N), by = method]
 test[, .(
   list(FGR(Hist(time, D) ~ X + Z, data = imp_dats, cause = 1))
 ), by = c("method", "imp_ind")]

 microbenchmark::microbenchmark(
   "lapply" = {
     purrr::modify_depth(imp_dats, .depth = 1L, .f = ~ {
       lapply(.x, function(imp) FGR(model_formula, data = imp, cause = 1))
     })
   },
   "nested" = {
     nested_imps_mods <- nested_imps[, .(
       mods = list(
         lapply(imps[[1]], function(imp_dat) FGR(model_formula, data = imp_dat, cause = 1))
       )
     ), by = method]
   }, times = 3
 )

  # Bind the true ones..too
  # intersect() colnames for binding
}
