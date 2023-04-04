source(here("packages.R"))
source(here("R/data-generation.R"))

df <- generate_dataset(
  n = 1000,
  args_event_times = list(
    mechanism = "correct_FG",
    params = tar_read(true_params_correct_FG),
    censoring_type = "exponential",
    censoring_params = list("exponential" = 0.46)
  ),
  args_missingness = list(mech_params = list("prob_missing" = 0.1, "mechanism_expr" = "Z"))
)


# Prep args
originaldata <- data.frame(df)
smtype <- "coxph"
smformula <- "Surv(time, D) ~ X + Z"
method <- make.method(originaldata, defaultMethod = c("norm", "logreg", "mlogit", "podds"))
m <- 3
numit <- 2
cause <- 1
kmi_args <- list(
  formula = Surv(time, D != 0) ~ 1,
  data = originaldata,
  etype = "D",
  failcode = cause,
  nimp = m
)

# Add parameters for MI imps per kmi imp??
smcfcs.finegray <- function(originaldata,
                            smtype = "coxph",
                            smformula,
                            method,
                            m = 5,
                            numit = 10,
                            rjlimit = 5000,
                            #kmi_args = NULL,
                            ...) {

  # Or maybe avoid this and force D to be numeric..
  # CHeck also how kmi does this part
  outcome_vars <- all.vars(update(as.formula(smformula), . ~ 1))
  status_var_name <- outcome_vars[length(outcome_vars)]
  status_var <- originaldata[, status_var_name]
  cens_code <- ifelse(is.factor(status_var), levels(status_var)[1], 0L)

  form_split <- unlist(strsplit(x = smformula, split = "~"))
  surv_obj <- eval(parse(text = form_split[1]), envir = originaldata |>
                     mutate(D = as.numeric(D)))
  attr()

  form <- as.formula(smformula)
  # use crossprod in data generation?
  x <- eval(
    update(as.formula(smformula), . ~ 0),
    envir = originaldata
  )

  #debugonce(kmi)
  kmi_imps <- kmi::kmi(
    Surv(time, D != 0) ~ 1, # or > 0, add cens code
    data = data.frame(originaldata),
    etype = substitute(D),
    failcode = 1,
    nimp = m
  )

  kmi_imps

  # Do invisible and capture output bits here..

}
