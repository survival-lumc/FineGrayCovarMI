# Add parameters for number of smcfcs imps per kmi imp?
smcfcs.finegray <- function(originaldata,
                            smformula,
                            method,
                            cause = 1,
                            m = 5,
                            numit = 10,
                            rjlimit = 5000,
                            kmi_args = list(
                              "formula" = ~ 1,
                              "bootstrap" = FALSE,
                              "nboot" = 10,
                              "epsilon" = 1
                            ),
                            ...) {

  # Locate/sort out outcome variables
  outcome_vars <- all.vars(update(as.formula(smformula), . ~ 1))
  time_var_name <- head(outcome_vars, n = 1L) # will need to add warning/error if (tstart, tstop) format
  status_var_name <- tail(outcome_vars, n = 1L)
  time_var <- originaldata[[time_var_name]]
  status_var <- originaldata[[status_var_name]]
  if (!is.numeric(status_var))
    stop("Status variable should be numeric, with 0 indicating a censored observation.")

  # Get functional form RHS of formula
  smformula_rhs <- unlist(strsplit(smformula, split = "~"))[2]

  # Use names consistent with kmi's names
  smformula_processed <- paste0("Surv(newtimes, newevent) ~", smformula_rhs)
  meths_smcfcs <- c(method, "newtimes" = "", newevent = "")

  # Check number censored, since based on this we use {kmi} or not
  num_censored <- sum(status_var == 0)

  # If no censoring: just pre-process data
  if (num_censored == 0) {

    # Set competing events to event time larger that largest event 1 time
    # (in this case: largest event 1 time + 0.1)
    eps <- 0.1
    max_ev1_time <- max(time_var[status_var == cause])
    newtimes <- time_var
    newtimes[status_var != cause] <- max_ev1_time + eps
    newevent <- as.numeric(status_var == cause)

    # Bind to original data
    originaldata_processed <- cbind.data.frame(originaldata, newtimes, newevent)

    # Run imputations
    smcfcs_obj <- smcfcs(
      originaldata = originaldata_processed,
      smtype = "coxph",
      smformula = smformula_processed,
      method = meths_smcfcs,
      rjlimit = rjlimit,
      numit = numit,
      m = m,
      ...
    )

  } else { # Need to multiply impute censoring times!

    # Prepare kmi() formula (for now default Kaplan-Meier imputation)
    lhs_kmi <- paste0("Surv(", paste(outcome_vars, collapse = ", "), " != 0)")
    form_kmi <- as.formula(paste0(lhs_kmi, deparse1(kmi_args$formula)))
    args_cens_imps <- c(
      list(
        "formula" = form_kmi,
        "data" = originaldata,
        "etype" = as.symbol(status_var_name),
        "failcode" = cause,
        "nimp" = m
      ),
      kmi_args[-1] # remove formula
    )

    # Impute missing censoring times in first loop
    # We use the local timefixed version (in smcfcs package: edit survfit as globally new fun with timefix = FALSE)
    kmi_imps <- do.call(kmi, args = args_cens_imps)

    # And now impute the covariates
    imps_loop <- lapply(kmi_imps$imputed.data, function(new_outcomes) {

      df_imp <- cbind(kmi_imps$original.data, new_outcomes)
      df_imp$newevent <- as.numeric(df_imp$newevent) - 1L # Make numeric for smcfcs

      smcfcs_modif <- smcfcs(
        originaldata = df_imp,
        smtype = "coxph",
        smformula = smformula_processed,
        method = meths_smcfcs,
        rjlimit = rjlimit,
        numit = numit,
        m = 1L, # one imputation per kmi dataset
        ...
      )

      return(smcfcs_modif)
    })

    smcfcs_obj <- smcfcs:::combine_smcfcs_objects(imps_loop)
  }

  # Do invisible and capture output bits here?
  return(smcfcs_obj)
}
