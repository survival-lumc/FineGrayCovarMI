# New tests
library(data.table)
library(riskRegression)
library(survival)
library(mice)
library(prodlim)
library(kmi)
library(smcfcs)
library(tidyverse)
library(broom)

# set.seed(...)

source("R/data-generation.R")

params_direct <- list(
  "cause1" = list(
    "betas" = c(0.5, 0.25),
    "p" = 0.4,
    "base_rate" = 0.5,
    "base_shape" = 1
  ),
  "cause2" = list(
    "betas" = c(0.5, 0.25),
    "base_rate" = 0.15,
    "base_shape" = 1
  )
)

dat <- generate_complete_dataset(
  n = 2000,
  params = params_direct,
  model_type = "direct",
  X_type = "binary",
  predictor_formulas = list("cause1" = ~  X + Z, "cause2" = ~ X + Z),
  control = list("admin_cens_time" = 1e20, "cens_rate" = 1e-10)
)

dat[]
table(dat$D)

FGR(Hist(time, D) ~ X + Z, data = dat, cause = 1)
coxph(Surv(time, D == 2) ~ X + Z, data = dat)


# Brief sims

map_dfr(
  .x = seq_len(250L),
  .f = ~ {
    dat <- generate_complete_dataset(
      n = 1000,
      params = params_direct,
      model_type = "direct",
      X_type = "continuous", # binary
      predictor_formulas = list("cause1" = ~ X + Z, "cause2" = ~ X + Z),
      control = list("admin_cens_time" = 1e20, "cens_rate" = 1e-10) # add option of NO censoring
    )
    mod_fg <- FGR(Hist(time, D) ~ X + Z, data = dat, cause = 1)$crrFit
    tidy(mod_fg)
    #mod_squeeze <- coxph(Surv(time, D == 2) ~ X + Z, data = dat)
    #bind_rows(tidy(mod_fg), tidy(mod_squeeze), .id = "mod_cause")
  },
  .id = "rep"
) |>
  #mutate(mod_cause = factor(mod_cause, labels = c("FG", "squeeze"))) |>
  ggplot(aes(term, estimate)) +
  geom_boxplot() +
  geom_hline(yintercept = c(0.5, 0.25), linetype = "dotted") #+
  #facet_wrap(~ mod_cause)


#
dat <- generate_complete_dataset(
  n = 50000,
  params = params_direct,
  model_type = "direct",
  X_type = "continuous", # binary
  predictor_formulas = list("cause1" = ~ X + Z, "cause2" = ~ X + Z),
  control = list("admin_cens_time" = 1e20, "cens_rate" = 1e-5) # add option of NO censoring
)
mod_fg <- FGR(Hist(time, D) ~ X + Z, data = dat, cause = 1)$crrFit
tidy(mod_fg)

#other fast packages
library(fastcmprsk)
mod_fg_fast <- fastCrr(Crisk(time, D, cencode = 0) ~ X + Z, variance = TRUE, data = dat)
mod_fg_fast$coef



# Simple MI no censoring --------------------------------------------------

set.seed(984984)
# Perhaps msprep quicker? but FGR much better for prediction

dat <- generate_complete_dataset(
  n = 2000,
  params = params_direct,
  model_type = "direct",
  X_type = "continuous", # binary
  predictor_formulas = list("cause1" = ~ X + Z, "cause2" = ~ X + Z),
  control = list("admin_cens_time" = 1e20, "cens_rate" = 1e-10) # add option of NO censoring
)
table(dat$D)

FGR(Hist(time, D) ~ X + Z, data = dat, cause = 1)
dat[, ':=' (
  D_star = as.numeric(D == 1),
  time_star = ifelse(D == 1, time, max(time) + 1e-8), # works better than setting to super large number..
  miss_ind = rbinom(.N, size = 1, prob = 0.4),
  X_compl = X
)]
coxph(Surv(time_star, D_star) ~ X_compl + Z, data = dat)
dat[, X := ifelse(miss_ind == 1, NA_real_, X)]
make.method(dat, defaultMethod = c("norm", "logreg", "mlogit", "podds"))
imps <- smcfcs(
  originaldata = data.frame(dat),
  smtype = "coxph",
  smformula = "Surv(time_star, D_star) ~ X + Z",
  method = make.method(dat, defaultMethod = c("norm", "logreg", "mlogit", "podds")),
  rjlimit = 5000,
  numit = 30,
  m = 10
)
plot(imps)
lapply(imps$impDatasets, function(imp) coxph(Surv(time_star, D_star) ~ X + Z, data = imp)) |>
  pool() |>
  tidy()

# With censoring ----------------------------------------------------------

set.seed(84645)

dat_cens <- generate_complete_dataset(
  n = 2000,
  params = params_direct,
  model_type = "direct",
  X_type = "continuous", # binary
  predictor_formulas = list("cause1" = ~ X + Z, "cause2" = ~ X + Z),
  control = list("admin_cens_time" = 1e20, "cens_rate" = 0.05) # add option of NO censoring
)

# Test KMI
table(dat_cens$D)
dat_cens[]
kmi_single <- kmi(
  Surv(time, D > 0) ~ 1,
  data = data.frame(dat_cens),
  etype = D,
  failcode = 1,
  nimp = 1
)
# The vector to impute
kmi_single$imputed.data[[1]][kmi_single$original.data$D == 2, ]$newtimes

# Try adjusting
dat_cens[, ':=' (
  t_star = ifelse(D == 2, NA_real_, time)#,
  #D_star = as.numeric(D == 1)
)]
dat_cens2 <- copy(dat_cens)
#dat_cens2[, time := NULL]
dat_cens2

meths <- make.method(dat_cens2)
predmat <- make.predictorMatrix(dat_cens2)
predmat[] <- 0
predmat["t_star", c("time", "D")] <- 1
predmat
#forms <- make.formulas(
#  dat_cens2,

#)
#forms$t_star <- t_star ~ Z

impute_cens_times <- function(time, D) {

  browser()
  kmi_single <- kmi(
    Surv(time, D > 0) ~ 1,
    data = cbind.data.frame("time" = time, "D" = D),
    etype = D,
    failcode = 1,
    nimp = 1
  )
  # The vector to impute
  kmi_single$imputed.data[[1]][kmi_single$original.data$D == 2, ]$newtimes
}

#impute_cens_times(dat_cens2$time, dat_cens2$D)


meths["t_star"] <- paste(
  "~I(", expression(impute_cens_times(time, D)),")"
)
meths
imps <- mice(
  data = data.frame(dat_cens2),
  method = meths,
  m = 2,
  maxit = 2,
  predictorMatrix = predmat
)

imp_dat <- cbind(
  kmi_single$original.data,
  kmi_single$imputed.data[[1]]
)
imp_dat[imp_dat$D == 2, ] |> head()
imp_dat |>  head()

# Generate single KMI

table(dat_cens$D)
# Testing KMI
kmi_obj <- kmi(
  Surv(time, D > 0) ~ 1,
  data = data.frame(dat_cens),
  etype = D,
  failcode = 1,
  nimp = 2
)
kmi_obj

kmi_obj$imputed.data
kmi_obj$imputed.data[[1]]
