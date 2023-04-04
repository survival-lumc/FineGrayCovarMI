library(data.table)
library(riskRegression)
library(survival)
library(mice)
library(prodlim)
library(kmi)
library(smcfcs)
library(tidyverse)
library(mstate)
library(broom)

source("R/data-generation.R")

one_replication <- function() {

  # Parameters
  params_direct <- list(
    "cause1" = list(
      "betas" = c(0.75, 0.5),
      "p" = 0.4,
      "base_rate" = 0.5,
      "base_shape" = 1
    ),
    "cause2" = list(
      "betas" = c(0.75, 0.5),
      "base_rate" = 0.5,
      "base_shape" = 1
    )
  )

  # Global imp settings
  m <- 10 #15
  iters_smcfcs <- 20
  rjlimit_smcfcs <- 5000L

  # Generate complete data
  dat <- generate_complete_dataset(
    n = 1000,
    params = params_direct,
    model_type = "direct",
    X_type = "binary",
    predictor_formulas = list("cause1" = ~ X + Z, "cause2" = ~ X + Z),
    control = list("admin_cens_time" = NULL, "cens_rate" = 0.15) # 20% of obs are censored
  )

  # Add missings; note subdist and cause-spec diverge as time goes on
  process_pre_imputing(dat)

  # -- Start
  form <- Hist(time, D) ~ X + Z

  # Method 0: Full data
  mod_full <- FGR(Hist(time, D) ~ X_compl + Z, data = dat, cause = 1)

  # Method 1: CCA
  mod_CCA <- FGR(form, data = dat, cause = 1)

  # Method 2: MICE as if comp risks
  meths <- make.method(dat, defaultMethod = c("norm", "logreg", "polyreg", "polr"))
  predmat <- make.predictorMatrix(dat)
  predmat[] <- 0
  predmat["X", c("Z", "D", "H_cause1", "H_cause2")] <- 1

  mice_comp <- mice(
    data = data.frame(dat),
    method = meths,
    m = m,
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
    rjlimit = rjlimit_smcfcs,
    numit = iters_smcfcs,
    m = m
  )

  # Method 5: SMC-FCS w/ modified dataset (do in parallel later!)

  # Start with kmi step
  kmi_imps <- kmi(
    Surv(time, D != 0) ~ 1,
    data = data.frame(dat),
    etype = D,
    failcode = 1,
    nimp = m
  )

  # Intermediate methods objects
  meths_smcfcs <- make.method(
    cbind(dat, kmi_imps$imputed.data[[1]]),
    defaultMethod = c("norm", "logreg", "mlogit", "podds")
  )

  # Single smcfcs imp per kmi imputed dataset
  imp_dats_smcfcs <- lapply(kmi_imps$imputed.data, function(new_outcomes) {

    df_imp <- cbind(kmi_imps$original.data, new_outcomes)
    df_imp$newevent <- as.numeric(df_imp$newevent) - 1L # Just for smcfcs

    # Use switch for imp methods
    smcfcs_modif <- smcfcs(
      originaldata = df_imp,
      smtype = "coxph",
      smformula = "Surv(newtimes, newevent) ~ X + Z",
      method = meths_smcfcs,
      rjlimit = rjlimit_smcfcs,
      numit = iters_smcfcs,
      m = 1L
    )

    return(smcfcs_modif)
  })

  # Pool already here
  mods_smcfcs <- lapply(
    imp_dats_smcfcs,
    function(imp) coxph(Surv(newtimes, newevent) ~ X + Z, data = imp)
  )

  # Method 6: MICE w/ cumulative subdistribution hazard - (no actual need for kmi)
  predmat_subdist <- predmat
  predmat_subdist[] <- 0
  predmat_subdist["X", c("Z", "D_star", "H_subdist_cause1")] <- 1

  mice_subdist <- mice(
    data = data.frame(dat),
    method = meths,
    m = m,
    maxit = 1,
    predictorMatrix = predmat_subdist,
    printFlag = FALSE
  )

  # Store all imputed datasets
  imps_dats <- list(
    "mice_comp" = complete(mice_comp, action = "all"),
    "smcfcs_comp" = smcfcs_comp$impDatasets,
    #"smcfcs_modif" = smcfcs_modif$impDatasets,
    "mice_subdist" = complete(mice_subdist, action = "all")
  )

  rm(mice_comp, smcfcs_comp, mice_subdist, imp_dats_smcfcs) # probably will need to keep stuff for cuminc prediction..

  pooled_estims <- modify_depth(imps_dats, .depth = 1L, .f = ~ {
    mods <- lapply(.x, function(imp) FGR(form, data = imp, cause = 1)$crrFit)
    tidy(pool(mods))
  })

  res <- bind_rows(
    c(
      list("full" = tidy(mod_full$crrFit)),
      list("CCA" = tidy(mod_CCA$crrFit)),
      pooled_estims,
      list("SMC-FCS + kmi" = tidy(pool(mods_smcfcs)))
    ), .id = "method"
  ) |>
    mutate(term = fct_collapse(term, "X" = c("X_compl1", "X2"), "Z" = "Z"))
  return(res)
}

set.seed(8762469)
system.time({
  test <- map_dfr(seq_len(10L), .f = ~ {
    cat(paste0("\n\nReplication #", .x, "\n\n"))
    one_replication()
  })
})

# test |>
#   group_by(term, method) |>
#   summarise(
#     avg = mean(estimate),
#     emp_se = sd(estimate)
#   )

saveRDS(test, "sims_cens.rds")
