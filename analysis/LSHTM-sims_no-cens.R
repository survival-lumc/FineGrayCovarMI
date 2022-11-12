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

# set.seed(...)

source("R/data-generation.R")

# Parameters
params_direct <- list(
  "cause1" = list(
    "betas" = c(0.5, 0.25),
    "p" = 0.4,
    "base_rate" = 0.5,
    "base_shape" = 1
  ),
  "cause2" = list(
    "betas" = c(0.5, 0.25),
    "base_rate" = 0.5,
    "base_shape" = 1
  )
)

# Global imp settings
m <- 15
iters <- 20

# Go through one rep
dat <- generate_complete_dataset(
  n = 2000, # up to 2000
  params = params_direct,
  model_type = "direct",
  X_type = "binary", # or "continous"
  predictor_formulas = list("cause1" = ~ X + Z, "cause2" = ~ X + Z),
  control = list("admin_cens_time" = 1e20, "cens_rate" = 1e-10) # add option of NO censoring
)

table(dat$D)

# Add missing data
dat[, ':=' (
  time_star = ifelse(D == 1, time, max(time) + 1e-8), # works better than setting to super large number.. @Hein?
  D_star = as.numeric(D == 1),
  miss_ind = rbinom(.N, size = 1, prob = plogis(-0.5 + Z)),
  X_compl = X,
  H_cause1 = add_marginal_cumhaz(
    timevar = time,
    statusvar = D,
    cause = 1,
    type = "cause_spec"
  ),
  H_cause2 = add_marginal_cumhaz(
    timevar = time,
    statusvar = D,
    cause = 2,
    type = "cause_spec"
  ),
  H_subdist_cause1 = add_marginal_cumhaz(
    timevar = time,
    statusvar = D,
    cause = 1,
    type = "subdist"
  )
)]
dat[, ':=' (
  X = factor(ifelse(miss_ind == 1, NA_character_, X)), # NA_real_ if continuous
  D = factor(D), # for imp model
  H_modif_cause1 = add_marginal_cumhaz(
    timevar = time_star,
    statusvar = D_star,
    cause = 1,
    type = "cause_spec"
  )
)]
setorder(dat, "time")
dat[]
mean(dat$miss_ind)

# Change this to ggplot, against time (but not for modif one)
plot(cmprsk::cuminc(dat$time, dat$D))
#

# Method 1: CCA
form <- Hist(time, D) ~ X + Z
mod_CCA <- FGR(form, data = dat, cause = 1)
#tidy(mod_CCA$crrFit)

# Method 2: MICE as if comp risks
meths <- make.method(dat, defaultMethod = c("norm", "logreg", "polyreg", "polr"))
predmat <- make.predictorMatrix(dat)
predmat[] <- 0
predmat["X", c("Z", "D", "H_cause1", "H_cause2")] <- 1
meths
predmat

mice_comp <- mice(
  data = data.frame(dat),
  method = meths,
  m = m,
  maxit = iters,
  predictorMatrix = predmat
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
  rjlimit = 5000,
  numit = iters,
  m = m
)

# Method 4: MICE w/ modified dataset
predmat_modif <- predmat
predmat_modif[] <- 0
predmat_modif["X", c("Z", "D_star", "H_modif_cause1")] <- 1
predmat_modif

mice_modif <- mice(
  data = data.frame(dat),
  method = meths,
  m = m,
  maxit = iters,
  predictorMatrix = predmat_modif
)

# Method 5: SMC-FCS w/ modified dataset
smcfcs_modif <- smcfcs(
  originaldata = data.frame(dat),
  smtype = "coxph",
  smformula = "Surv(time_star, D_star) ~ X + Z",
  method = meths_smcfcs,
  rjlimit = 5000,
  numit = iters,
  m = m
)

# Method 6: MICE w/ cumulative subdistribution hazard
predmat_subdist <- predmat
predmat_subdist[] <- 0
predmat_subdist["X", c("Z", "D_star", "H_subdist_cause1")] <- 1
predmat_subdist

mice_subdist <- mice(
  data = data.frame(dat),
  method = meths,
  m = m,
  maxit = iters,
  predictorMatrix = predmat_subdist
)

# Check convergence of all
plot(mice_subdist)
plot(mice_comp)
plot(mice_modif)
plot(smcfcs_comp)
plot(smcfcs_modif)

# Store all imputed datasets
imps_dats <- list(
  "mice_comp" = complete(mice_comp, action = "all"),
  "smcfcs_comp" = smcfcs_comp$impDatasets,
  "mice_modif" = complete(mice_modif, action = "all"),
  "smcfcs_modif" = smcfcs_modif$impDatasets,
  "mice_subdist" = complete(mice_subdist, action = "all")
)

pooled_estims <- modify_depth(imps_dats, .depth = 1L, .f = ~ {
  mods <- lapply(.x, function(imp) FGR(form, data = imp, cause = 1)$crrFit)
  tidy(pool(mods))
})

res <- bind_rows(c(list("CCA" = tidy(mod_CCA$crrFit)), pooled_estims), .id = "method")




# Numerical example -------------------------------------------------------


# Why we need to re-weight..


