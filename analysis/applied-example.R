library(data.table)
library(survival)
library(smcfcs)
library(mice)
library(broom)
library(prodlim)
library(riskRegression)
library(future)
library(future.apply)
options(contrasts = rep("contr.treatment", 2))

#
set.seed(89648)


dat <- readRDS("data/dat_clean.rds")
setDT(dat)

outcome_vars <- c(
  "time_ci_adm",
  "status_ci_adm"
)

predictors <- c(
  "hctci_risk",
  "age_allo1_decades",
  "wbc_allo1",
  "hb_allo1",
  "pb_allo1",
  "sweat_allo1",
  "WEIGLOSS_allo1",
  "KARNOFSK_threecat",
  "PATSEX",
  "donrel_bin",
  "cmv_match",
  "ruxo_preallo1",
  "ric_allo1"
)

censoring <- c(
  "year_allo1",
  "intdiagtr_allo1",
  "tbi_allo1"
)

dat_sub <- dat[, c(
  "time_ci_adm",
  "status_ci_adm",
  "age_allo1_decades",
  "PATSEX",
  "hctci_risk", # make them ordered later!!
  "year_allo1",
  "intdiagtr_allo1",
  "KARNOFSK_threecat",
  "ruxo_preallo1",
  "donrel_bin",
  "cmv_match",
  "wbc_allo1", # what analysis model should we use??
  "ric_allo1"
)]
dat_sub <- dat_sub[complete.cases(status_ci_adm)]
dat_sub[, ':=' (
  KARNOFSK_threecat = as.ordered(KARNOFSK_threecat),
  hctci_risk = as.ordered(hctci_risk)
)]

naniar::gg_miss_upset(dat_sub)
naniar::miss_var_summary(dat_sub)
naniar::prop_complete_case(dat_sub)

# Make model for censoring!
mod_cens <- coxph(
  Surv(time_ci_adm, status_ci_adm == 0) ~ .,
  data = dat_sub
)
summary(mod_cens)

sm_form <- reformulate(
  response = "Surv(time_ci_adm, status_ci_adm)",
  termlabels = c(
    "age_allo1_decades",
    "PATSEX",
    "hctci_risk", # make them ordered later!!
    "year_allo1",
    "intdiagtr_allo1",
    "KARNOFSK_threecat",
    "donrel_bin",
    "ruxo_preallo1",
    "cmv_match",
    "wbc_allo1",
    "ric_allo1"
  )
)
meths_smcfcs <- mice::make.method(dat_sub, defaultMethod = c("norm", "logreg", "mlogit", "podds"))


source("R/smcfcs.finegray_edit.R")
source("R/kmi-timefixed.R")

imps <- smcfcs.finegray.edit(
  originaldata = data.frame(dat_sub),
  smformula = deparse1(sm_form),
  method = meths_smcfcs,
  cause = 2, #nrm
  m = 10,
  #numit = 20,
  #rjlimit = 10000,
  #kmi_args = "1" # just the rhs
  kmi_args = c(
    "1"
    #"year_allo1",
    #"intdiagtr_allo1"
  ) # just the rhs
)

plot(imps, scales = "free")

model_formula <- update(sm_form, Surv(newtimes, newevent) ~ .)

imp_mods <- lapply(
  imps$impDatasets,
  function(imp_dat) coxph(model_formula, data = imp_dat)
)

tidy(pool(imp_mods), exponentiate = TRUE)[, c("term", "estimate", "std.error", "p.value")]


# Try with complete cases:
imps_cca <- lapply(
  imps$impDatasets,
  function(imp_dat) {
    df <- merge(
      imp_dat[, c("newtimes", "newevent", "age_allo1_decades", "intdiagtr_allo1")],
      dat_sub
    )
    coxph(model_formula, data = df)
  }
)

naniar::prop_complete_case(dat_sub)

tidy(pool(imps_cca), exponentiate = TRUE)[, c("term", "estimate", "std.error", "p.value")]

imps_csh <- smcfcs(
  originaldata = data.frame(dat_sub),
  smformula = list(
    deparse1(update(sm_form, Surv(time_ci_adm, status_ci_adm == 1) ~ .)),
    deparse1(update(sm_form, Surv(time_ci_adm, status_ci_adm == 2) ~ .))
  ),
  m = 10,
  method = meths_smcfcs,
  smtype = "compet",
  rjlimit = 5000
)

plot(imps_csh, scales = "free")

# - Checking prop assumption?
# - Also what do we plot against each other?
# - Check with the admin censoring what is up?

predt <- 24 # months

preds_mi_fg <- lapply(
  imps$impDatasets,
  function(imp) {
    mod <- do.call(coxph, args = list(formula = model_formula, data = imp, x = TRUE))
    data.table(
      pred = drop(predictRisk(mod, newdata = imp, times = predt)),
      dat_sub
    )
  }
)

preds_mi_cs <- lapply(
  imps_csh$impDatasets,
  function(imp) {
    mod <- do.call(
      riskRegression::FGR,
      args = list(
        formula = update(model_formula, Hist(time_ci_adm, status_ci_adm) ~ .),
        data = imp,
        cause = 2
      )
    )
    data.table(
      pred = drop(predictRisk(mod, newdata = imp, times = predt)),
      dat_sub
    )
  }
)

# Plot predictions against each outer??
predos <- rbind(
  cbind(rbindlist(preds_mi_fg), "imp_method" = "smcfcs_fg"),
  cbind(rbindlist(preds_mi_cs), "imp_method" = "smcfcs_csh")
)

predos <- cbind(rbindlist(preds_mi_fg), "pred_csh" = rbindlist(preds_mi_cs)$pred)

# Hmm hard to compare because imputed values must vary so much??
hize <- predos[, .(
  .N,
  fg = mean(pred),
  pred_pooled_csh = mean(pred_csh)
), by = c(
  "time_ci_adm",
  "age_allo1_decades",
  "intdiagtr_allo1"
)]

plot(hize$fg, hize$pred_pooled_csh, xlim = c(0,1), ylim = c(0, 1))


library(ggplot2)

# ADD ID AGAIN!!
hi <- predos[, .(
  .N,
  pred_pooled = mean(pred)
), by = c(
  "imp_method",
  "time_ci_adm",
  "age_allo1_decades",
  "intdiagtr_allo1"
)]

plot(
  hi[imp_method == "smcfcs_fg"][["pred_pooled"]],
  hi[imp_method == "smcfcs_csh"][["pred_pooled"]],
  xlim = c(0, 1),
  ylim = c(0, 1)
)

hi |>
  ggplot(aes(pred_pooled, fill = imp_method)) +
  geom_density()
