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

set.seed(123449)

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
m <- 50 #15
iters_smcfcs <- 20
rjlimit_smcfcs <- 5000L

# Generate complete data
dat <- generate_complete_dataset(
  n = 2000, # up to 2000
  params = params_direct,
  model_type = "direct",
  X_type = "binary", # or "continous"
  predictor_formulas = list("cause1" = ~ X + Z, "cause2" = ~ X + Z),
  control = list("admin_cens_time" = 1e20, "cens_rate" = 1e-10) # add option of NO censoring
)

# Add censoring...
dat[, cens := rexp(.N, rate = 0.1)]
dat[, ':=' (
  time = pmin(time, cens),
  D = ifelse(time < cens, D, 0L)
)]
table(dat$D)
modo <- FGR(Hist(time, D) ~ X + Z, data = dat, cause = 1)
modo$crrFit$coef

dat[,  miss_ind := rbinom(.N, size = 1, prob = plogis(-0.5 + Z))]
dat[, X := factor(ifelse(miss_ind == 1, NA_character_, X))]
dat

# Run kmi imps
kmi_imps <- kmi(
  Surv(time, D > 0) ~ 1,
  data = data.frame(dat),
  etype = D,
  failcode = 1,
  nimp = m
)

cox.kmi(
  Surv(time, D == 1) ~ X  + Z,
  kmi_imps
) |> summary()

kmi_imps$imputed.data[[1]]$newevent

meths_smcfcs <- make.method(
  cbind(dat, kmi_imps$imputed.data[[1]]),
  defaultMethod = c("norm", "logreg", "mlogit", "podds")
)

imp_dats <- lapply(kmi_imps$imputed.data, function(new_outcomes) {

  df_imp <- cbind(kmi_imps$original.data, new_outcomes)
  df_imp$newevent <- as.numeric(df_imp$newevent) - 1L

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

  return(smcfcs_modif$impDatasets[[1]])
})

lapply(
  imp_dats,
  function(imp) coxph(Surv(newtimes, newevent) ~ X + Z, data = imp)
) |>
  pool() |>
  tidy()

# Add missings
#dat_to_impute <- process_pre_imputing(dat)

# Settings to put in function: number of kmi imps, number of smcfcs *per* kmi imp
