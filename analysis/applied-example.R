library(data.table)
library(survival)
library(smcfcs)
library(mice)
library(future)
library(future.apply)
library(riskRegression) # for predictRisk()
library(kmi)
library(ggplot2)
library(patchwork)

# Global seed
options(contrasts = rep("contr.treatment", 2))
set.seed(123978)

# Source-in necessary functions
source("R/data-generation.R")


# Data set-up -------------------------------------------------------------


# Directly from Polverelli imputation set-up
# Dataset created on 20/07/2023, cohort 2009-2019
dat <- data.table(readRDS("data/dat_clean.rds"))

# Set-up predictors and auxiliary
sm_predictors <- c(
  "age_allo1_decades", # this one is centered
  "PATSEX",
  "hctci_risk",
  "KARNOFSK_threecat",
  "ruxo_preallo1",
  "donrel_bin",
  "cmv_match",
  "wbc_allo1",
  "hb_allo1",
  "ric_allo1",
  "WEIGLOSS_allo1"
)

aux_predictors <- c(
  "year_allo1",
  "intdiagtr_allo1",
  "tbi_allo1"
)

# Subset of relevant vars
dat_sub <- dat[, c(
  # Outcomes
  "time_ci_adm",
  "status_ci_adm",
  sm_predictors,
  aux_predictors
), with = FALSE]

# Other possible disease vars:
# pb_allo1, sweat_allo1
# Include weight loss, but avoid pb since percentage..

# Use complete outcome data, and add id column
dat_sub <- dat_sub[complete.cases(status_ci_adm), ]
dat_sub[, id := seq_len(.N)]

# Edit some variables pre-imputation
dat_sub[, ':=' (
  KARNOFSK_threecat = as.ordered(KARNOFSK_threecat),
  hctci_risk = as.ordered(hctci_risk),
  wbc_allo1 = log(wbc_allo1 + 1) # since very skewed
)]

# Some visuals and stats
naniar::gg_miss_upset(dat_sub)
naniar::miss_var_summary(dat_sub)
naniar::prop_complete_case(dat_sub)
table(dat_sub$status_ci_adm)

# For fun, scale age in decades per 100 years!
#hist(dat_sub$age_allo1_decades)
#dat_sub[, age_allo1_decades := age_allo1_decades / 10]
# Only matters for vars with missing vals

# On the censoring... -----------------------------------------------------


# Exploratory model for censoring, use mainly fully observed vars
mod_cens <- coxph(
  Surv(time_ci_adm, status_ci_adm == 0) ~ year_allo1 +
    I(intdiagtr_allo1 / 12) + # years
    donrel_bin +
    age_allo1_decades +
    PATSEX +
    cmv_match +
    ric_allo1 +
    KARNOFSK_threecat +
    tbi_allo1,
  data = dat_sub#[!(year_allo1 %in% c(8, 9, 10))] # Reduce to those with 5y potential follow-up
)

# TBI, RIC, Donor type should probs be included in model too..
summary(mod_cens)

# KM
sf_cens <- survfit(
  Surv(time_ci_adm, status_ci_adm != 0) ~ factor(year_allo1),
  data = dat_sub
)
ggsurvfit::ggsurvfit(sf_cens, linetype_aes = TRUE)


# First impute the censoring times for all methods ------------------------


# For estimation: quicker to impute cens times for all imp methods!


# Add the cumulative hazards first
dat_sub[, ':=' (
  H1 = compute_marginal_cumhaz(
    type = "cause_spec",
    timevar = time_ci_adm,
    statusvar = status_ci_adm,
    cause = 1
  ),
  H2 = compute_marginal_cumhaz(
    type = "cause_spec",
    timevar = time_ci_adm,
    statusvar = status_ci_adm,
    cause = 2
  ),
  # Add separate indicators
  D1 = as.numeric(status_ci_adm == 1),
  D2 = as.numeric(status_ci_adm == 2)
)]

# Check columns, and the hazards (would all cause hazard work well too here?)
colnames(dat_sub)

# Non-parametric cumincs, quite different to MDS example!
plot(cmprsk::cuminc(dat_sub$time_ci_adm, dat_sub$status_ci_adm))

# Cumulative hazards against each other
plot(dat_sub$time_ci_adm, dat_sub$H2)
points(dat_sub$time_ci_adm, dat_sub$H1)

# Factor for imputation of cens times
dat_sub$year_allo1 <- as.factor(dat_sub$year_allo1)

# Impute censoring times!!
cens_imps <- kmi(
  formula = Surv(time_ci_adm, status_ci_adm != 0) ~ year_allo1,
  data = data.frame(dat_sub),
  etype = status_ci_adm,
  failcode = 2, # non-relapse mortality is the outcome
  nimp = 10 # make bigger later, set globally
) # If we are running once, why not bootstrap?

# Checks with single imps, also to set-up imp methods
dat_single_imp <- cbind.data.frame(cens_imps$original.data, cens_imps$imputed.data[[1]])

# For set-up
dat_single_imp$H1_subdist <- compute_marginal_cumhaz(
  type = "cause_spec",
  timevar = dat_single_imp$newtimes,
  statusvar = dat_single_imp$newevent,
  cause = 2
)


# Set-up imputations ------------------------------------------------------


# First we need to know what the analysis model is
sm_form <- reformulate(
  response = "Surv(time_ci_adm, status_ci_adm)", termlabels = sort(sm_predictors)
)

# Set-up methods using just this first imp dataset
meths_smcfcs <- make.method(dat_single_imp, defaultMethod = c("norm", "logreg", "mlogit", "podds"))
meths_mice <- make.method(dat_single_imp) # pmm?

# Predictor matrices for mice
predmat_csh <- predmat_subdist <- make.predictorMatrix(dat_single_imp)
predmat_csh[, setdiff(colnames(predmat_csh), c(sm_predictors, "D1", "D2", "H1", "H2"))] <- 0
predmat_subdist[, setdiff(colnames(predmat_subdist), c(sm_predictors, "newevent", "H1_subdist"))] <- 0


# Now applying across imp datasets ----------------------------------------


# Set general imputation settings
n_imps <- 1 # because 1 imp per kmi dataset
n_cycles <- 10

# First store the imputed datasets together (and add subdist cumhaz)
cens_imp_dats <- lapply(
  X = cens_imps$imputed.data,
  FUN = function(imp_times) {
    dato <- cbind.data.frame(cens_imps$original.data, imp_times)
    dato$H1_subdist <- compute_marginal_cumhaz(
      type = "cause_spec",
      timevar = dato$newtimes,
      statusvar = dato$newevent,
      cause = 2
    )
    # Back to numeric
    dato$year_allo1 <- as.numeric(dat_sub$year_allo1)
    return(dato)
  }
)

# First actually check the variance in the subdist marginal haz
rbindlist(cens_imp_dats, idcol = "imp") |>
  ggplot(aes(newtimes, H1_subdist, group = factor(imp), col = factor(imp))) +
  geom_step() +
  theme(legend.position = "none")

# Hardly any difference! Just compute once?


# Run the covariate MIs.. -------------------------------------------------


# For parallel
rjlimit <- 10000
plan(multisession, workers = 3)

impdats_mice_csh <- future_lapply(
  X = cens_imp_dats,
  FUN = function(imp_times) {
    imps <- mice(
      data = imp_times,
      m = n_imps,
      maxit = n_cycles,
      method = meths_mice,
      predictorMatrix = predmat_csh,
      printFlag = FALSE
    )
    complete(imps, "all")[[1]]
  },
  future.seed = TRUE
)

impdats_mice_subdist <- future_lapply(
  X = cens_imp_dats,
  FUN = function(imp_times) {
    imps <- mice(
      data = imp_times,
      m = n_imps,
      maxit = n_cycles,
      method = meths_mice,
      predictorMatrix = predmat_subdist,
      printFlag = FALSE
    )
    complete(imps, "all")[[1]]
  },
  future.seed = TRUE
)

impdats_smcfcs_subdist <- future_lapply(
  X = cens_imp_dats,
  FUN = function(imp_times) {
    imp_times$newevent <- as.numeric(imp_times$newevent) - 1L
    imps <- smcfcs(
      originaldata = imp_times,
      smtype = "coxph",
      smformula = deparse1(update(sm_form, Surv(newtimes, newevent) ~ .)),
      m = n_imps,
      numit = n_cycles,
      method = meths_smcfcs,
      rjlimit = 5000
    )
    imp_dat <- imps$impDatasets[[1]]
    imp_dat$newevent <- factor(imp_dat$newevent, labels = c("0", "2"))
    imp_dat
  },
  future.seed = TRUE
)

impdats_smcfcs_csh <- future_lapply(
  X = cens_imp_dats,
  FUN = function(imp_times) {
    imps <- smcfcs(
      originaldata = imp_times,
      smtype = "compet",
      smformula = list(
        deparse1(update(sm_form, Surv(time_ci_adm, status_ci_adm == 1) ~ .)),
        deparse1(update(sm_form, Surv(time_ci_adm, status_ci_adm == 2) ~ .))
      ),
      m = n_imps,
      numit = n_cycles,
      method = meths_smcfcs,
      rjlimit = 5000
    )
    imp_dat <- imps$impDatasets[[1]]
    imp_dat$newevent <- as.numeric(imp_dat$newevent == 2)
    imp_dat
  },
  future.seed = TRUE
)

plan(sequential)



# Try stacking  -----------------------------------------------------------


imp_dats <- rbind(
  cbind(rbindlist(impdats_mice_subdist, idcol = "imp"), "method" = "MI subdist"),
  cbind(rbindlist(impdats_mice_csh, idcol = "imp"), "method" = "MICE CSH"),
  cbind(rbindlist(cens_imp_dats, idcol = "imp"), "method" = "CCA"),
  cbind(rbindlist(impdats_smcfcs_subdist, idcol = "imp"), "method" = "SMC-FCS subdist")
)

mods_imp_dats <- imp_dats[, .(
  mods = list(coxph(update(sm_form, Surv(newtimes, newevent == 2) ~ .), data = .SD))
), by = c("imp", "method")]

summ <- rbindlist(
  with(
    mods_imp_dats[, .(
      summary = list(broom::tidy(pool(mods)))#, exponentiate = TRUE, conf.int = TRUE))
    ), by = "method"],
    Map(f = cbind, method = method, summary)
  )
)

summ |>
  ggforestplot::forestplot(
    name = term,
    se = std.error,
    col = method,
    shape = method,
    logodds = TRUE,
    xtickbreaks = c(0.75, 1, 1.25, 1.5, 1.75),
    #xlim = c(0.5, 2.5),
    xlab = "Hazard ratio (95% CI)"
  )



# New predictions ---------------------------------------------------------


# Compare predictions on cases with complete data?
id_cc <- dat_sub[complete.cases(dat_sub), ][["id"]]
#
dat_sub[, prop_complete := 1 - rowSums(is.na(dat_sub)) / length(sm_predictors)]
horizon <- 24

preds_imp_dats <- imp_dats[, .(
  id = id,
  preds = {
    mod <- coxph(update(sm_form, Surv(newtimes, newevent == 2) ~ .), data = .SD, x = TRUE)
    drop(predictRisk(mod, newdata = .SD, times = horizon))
  }
), by = c("imp", "method")]

cloglog <- function(x) log(-log(1 - x))
inv_cloglog <- function(x) 1 - exp(-exp(x))

preds_cc <- preds_imp_dats[, .(
  preds_pooled = inv_cloglog(mean(cloglog(preds)))
), by = c("method", "id")][!(id %in% id_cc)]

wide <- dcast(preds_cc, id ~ method, value.var = "preds_pooled")
wide[, wt := dat_sub$prop_complete[match(id, dat_sub$id)]]

wide |>
  ggplot(aes(`MI subdist`, `SMC-FCS subdist`, col = wt)) +
  geom_point(alpha = 0.5, size = 3) +
  scale_color_viridis_b(direction = -1)

GGally::ggpairs(wide, columns = 3:5)

preds_cc |>
  ggplot(aes(preds_pooled, preds_pooled)) +
  geom_point() +
  facet_grid(~ method)



# Try a calibration ting
score_csh <- Score(
  list(test$preds_csh),
  data = cens_imp_dats[[1]], # because of ordering, use first cens imp dat
  times = horizon,
  formula = Hist(time_ci_adm, status_ci_adm) ~ 1,
  cause = 2,
  plots = "calibration"
)

#plotCalibration(score_csh, cens.method = "local", method = "quantile")
plotCalibration(score_csh, cens.method = "local")




# Try score using FGR/CSC on complete cases -------------------------------

library(prodlim)
mod_FGR <- FGR(
  update(sm_form, Hist(time_ci_adm, status_ci_adm) ~ .),
  data = dat_sub,
  cause = 2
)
score_cc <- Score(
  list(mod_FGR),
  data = dat_sub, # because of ordering, use first cens imp dat
  times = 2,
  formula = Hist(time_ci_adm, status_ci_adm) ~ 1,
  cause = 2,
  plots = "calibration"
)

plotCalibration(score_cc, cens.method = "local", method = "quantile")





