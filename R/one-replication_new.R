args_event_times <- list(
  mechanism = "correct_FG",
  params = tar_read(true_params_correct_FG),
  censoring_type = "none"
)
args_missingness <- list(mech_params = list("prob_missing" = 0.1, "mechanism_expr" = "Z"))
args_imputations <- list(m = 2, iters = 2, rjlimit = 1000)

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
  #...
}
