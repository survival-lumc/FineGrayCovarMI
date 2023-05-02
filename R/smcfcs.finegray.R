# Add parameters for MI imps per kmi imp??
smcfcs.finegray <- function(originaldata,
                            smformula,
                            method,
                            cause = 1,
                            m = 5,
                            numit = 10,
                            rjlimit = 5000,
                            #kmi_args = NULL,
                            ...) {

  # Locate/sort out outcome variables
  outcome_vars <- all.vars(update(as.formula(smformula), . ~ 1))
  time_var_name <- head(outcome_vars, n = 1L) # will need to add warning/error if tstart tstop
  status_var_name <- tail(outcome_vars, n = 1L)

  time_var <- originaldata[[time_var_name]]
  status_var <- originaldata[[status_var_name]]
  if (!is.numeric(status_var))
    stop("Status variable should be numeric, with 0 indicating a censored observation.")

  # Get function form RHS of formula
  smformula_rhs <- unlist(strsplit(smformula, split = "~"))[2]

  # Check number censored
  num_censored <- sum(status_var == 0)
  smformula_processed <- paste0("Surv(newtimes, newevent) ~", smformula_rhs)
  meths_smcfcs <- c(method, "newtimes" = "", newevent = "")

  # If no censoring: just pre-process data
  if (num_censored == 0) {

    # Naming here consistent with what kmi outputs
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
      m = m
    )

  } else {

    # Prepare kmi() formula (for now default kaplan-meier imp)
    lhs_kmi <- paste0("Surv(", paste(outcome_vars, collapse = ", "), " != 0)")
    form_kmi <- reformulate(termlabels = c("1"), response = lhs_kmi)

    # Add option to bootstrap?? See section 5.3.2 Beyersmann book

    # Impute missing censoring times in first loop
    kmi_imps <- do.call(
      kmi::kmi,
      args = list(
        "formula" = form_kmi,
        "data" = originaldata, # watch out here later
        "etype" = as.symbol(status_var_name),
        "failcode" = cause,
        "nimp" = m
      )
    )

    # Loop
    imps_loop <- lapply(kmi_imps$imputed.data, function(new_outcomes) {

      df_imp <- cbind(kmi_imps$original.data, new_outcomes)
      df_imp$newevent <- as.numeric(df_imp$newevent) - 1L # Just for smcfcs

      # Use switch for imp methods
      smcfcs_modif <- smcfcs(
        originaldata = df_imp,
        smtype = "coxph",
        smformula = smformula_processed,
        method = meths_smcfcs,
        rjlimit = rjlimit,
        numit = numit,
        m = 1L
      )

      return(smcfcs_modif)
    })

    smcfcs_obj <- smcfcs:::combine_smcfcs_objects(imps_loop)
  }

  # Do invisible and capture output bits here?
  return(smcfcs_obj)
}
