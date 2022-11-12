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
    n = 1000, # up to 2000
    params = params_direct,
    model_type = "direct",
    X_type = "binary", # or "continous"
    predictor_formulas = list("cause1" = ~ X + Z, "cause2" = ~ X + Z),
    control = list("admin_cens_time" = 1e20, "cens_rate" = 1e-10) # add option of NO censoring
  )

  # Add missings
  dat_to_impute <- process_pre_imputing(dat)

  # -- Start
  form <- Hist(time, D) ~ X + Z

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

  # Method 4: MICE w/ modified dataset
  predmat_modif <- predmat
  predmat_modif[] <- 0
  predmat_modif["X", c("Z", "D_star", "H_modif_cause1")] <- 1

  mice_modif <- mice(
    data = data.frame(dat),
    method = meths,
    m = m,
    maxit = 1,
    predictorMatrix = predmat_modif,
    printFlag = FALSE
  )

  # Method 5: SMC-FCS w/ modified dataset (do in parallel later!!!)
  smcfcs_modif <- smcfcs(
    originaldata = data.frame(dat),
    smtype = "coxph",
    smformula = "Surv(time_star, D_star) ~ X + Z",
    method = meths_smcfcs,
    rjlimit = rjlimit_smcfcs,
    numit = iters_smcfcs,
    m = m
  )

  # Method 6: MICE w/ cumulative subdistribution hazard
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
    "mice_modif" = complete(mice_modif, action = "all"),
    "smcfcs_modif" = smcfcs_modif$impDatasets,
    "mice_subdist" = complete(mice_subdist, action = "all")
  )

  rm(mice_comp, smcfcs_comp, mice_modif, smcfcs_modif, mice_subdist) # probably will need to keep stuff for cuminc prediction..

  pooled_estims <- modify_depth(imps_dats, .depth = 1L, .f = ~ {
    mods <- lapply(.x, function(imp) FGR(form, data = imp, cause = 1)$crrFit)
    tidy(pool(mods))
  })

  res <- bind_rows(c(list("CCA" = tidy(mod_CCA$crrFit)), pooled_estims), .id = "method")
  return(res)
}

system.time({
  test <- map_dfr(seq_len(2L), .f = ~ {
    cat(paste0("\n\nReplication #", .x, "\n\n"))
    one_replication()
  })
})

saveRDS(test, "sims_nocens.rds")

# test |>
#   ggplot(aes(term, estimate, fill = method)) +
#   geom_boxplot() +
#   theme_minimal() +
#   geom_hline(yintercept = c(0.5, 0.25))
#
# test |>
#   group_by(term, method) |>
#   summarise(avg = mean(estimate),
#             emp_se = sd(estimate))


# Just use crprep for everything innit;
