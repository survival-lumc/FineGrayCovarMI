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
#source("R/kmi-timefixed.R")

# Global seed
options(contrasts = rep("contr.treatment", 2))
#set.seed(89648)

source("R/data-generation.R")
source("R/crprep-timefixed.R")


# Data set-up -------------------------------------------------------------


# Directly from Polverelli imputation set-up
dat <- readRDS("data/dat_clean.rds")
setDT(dat)

# Subset of relevant vars
dat_sub <- dat[, c(
  # Outcomes
  "time_ci_adm",
  "status_ci_adm",
  # Covariates for predictions
  "age_allo1_decades", # this one is centered
  "PATSEX",
  "hctci_risk",
  "KARNOFSK_threecat",
  "ruxo_preallo1",
  "donrel_bin",
  "cmv_match",
  "wbc_allo1",
  "ric_allo1",
  # Good auxiliary variables
  "year_allo1",
  "intdiagtr_allo1",
  "tbi_allo1"
)]

# Other possible disease vars:
# wbc_allo1, hb_allo1, pb_allo1, sweat_allo1, WEIGLOSS_allo1..

# Use complete outcome data, and add id column
dat_sub <- dat_sub[complete.cases(status_ci_adm), ]
dat_sub[, id := seq_len(.N)]

# Set some variables as ordered for imputation
dat_sub[, ':=' (
  KARNOFSK_threecat = as.ordered(KARNOFSK_threecat),
  hctci_risk = as.ordered(hctci_risk)
)]

# Some visuals
naniar::gg_miss_upset(dat_sub)
naniar::miss_var_summary(dat_sub)
naniar::prop_complete_case(dat_sub)


# On the censoring... -----------------------------------------------------


# Dataset created on 20/07/2023, cohort 2009-2019

# Event props:
table(dat_sub$status_ci_adm)

# First start with the admin cens indicator
dat_sub[, ind_cens_adm := as.numeric(time_ci_adm == 60)]
table(dat_sub$ind_cens_adm) # 724 pats admin censored at 60 mo.s

# Pretend adm cens is 'cause 3' so that the censoring model is 'right'
dat_sub[, ind_cens_true := ifelse(ind_cens_adm == 1, 3, status_ci_adm)]

# 1292 of 'actual' censoring, including those Tx'ed in 2018-2019 which
# have not quite had 60 mo.s potential follow-up up
with(dat_sub, time_ci_adm[ind_cens_true == 0]) |> length()
dat_sub[year_allo1 == 10][order(time_ci_adm, decreasing = TRUE), .(
  time_ci_adm, status_ci_adm
)] #|> View()

# Exploratory model for 'true' censoring, use mainly fully observed vars
mod_cens <- coxph(
  Surv(time_ci_adm, ind_cens_true == 0) ~ year_allo1 +
    I(intdiagtr_allo1 / 12) + # years
    donrel_bin +
    age_allo1_decades +
    PATSEX +
    cmv_match +
    ric_allo1 +
    KARNOFSK_threecat +
    tbi_allo1,
  data = dat_sub[!(year_allo1 %in% c(8, 9, 10))] # Reduce to those with 5y potential follow-up
)

# TBI, RIC, Donor type should probs be included in model too..
summary(mod_cens)

# Note for KMI:
# - Difficult to specify mix of admin and random cens; if largest time
# is event time, it is added to support of cens distribution
# - Proposal: keep all admin + cens times together, and model to impute =
# continuous year of allo, continuous int diag, TBI, RIC, donor type
# - Not perfect, but should be granular enough
# - there 0.5-1% missing vals in TBI, RIC, donor type are dealt with
# internally by KMI using mean imp (proposal to accept and move on)


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
  #status_ci_adm = factor(status_ci_adm),
  # Add separate indicators
  D1 = as.numeric(status_ci_adm == 1),
  D2 = as.numeric(status_ci_adm == 2),
  wbc_allo1 = log(wbc_allo1 + 1) # because very skewed
)]

# Check columns, and the hazards (would all cause hazard work well too here?)
colnames(dat_sub)

# Non-parametric cumincs, quite different to MDS example!
plot(cmprsk::cuminc(dat_sub$time_ci_adm, dat_sub$status_ci_adm))

# Cumulative hazards against each other
plot(dat_sub$H1, dat_sub$H2)
abline(a = 0, b = 1)

plot(dat_sub$time_ci_adm, dat_sub$H2)
points(dat_sub$time_ci_adm, dat_sub$H1)

# Impute censoring times!!
cens_imps <- kmi(
  formula = Surv(time_ci_adm, status_ci_adm != 0) ~
    year_allo1 + intdiagtr_allo1 + ric_allo1 + tbi_allo1 + donrel_bin,
  data = data.frame(dat_sub),
  etype = status_ci_adm,
  failcode = 2, # non-relapse mortality is the outcome
  nimp = 25 # make bigger later
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
  response = "Surv(time_ci_adm, status_ci_adm)",
  termlabels = c(
    "age_allo1_decades",
    "PATSEX",
    "hctci_risk",
    "KARNOFSK_threecat",
    "donrel_bin",
    "ruxo_preallo1",
    "cmv_match",
    "wbc_allo1",
    "ric_allo1",
    "year_allo1", # include these two in analysis model or not? Data dictionary?
    "intdiagtr_allo1"
  )
)

sm_predictors <- all.vars(update(sm_form, 1 ~ .))

# Set-up methods using just this first imp dataset
meths_smcfcs <- make.method(dat_single_imp, defaultMethod = c("norm", "logreg", "mlogit", "podds"))
meths_mice <- make.method(dat_single_imp) # pmm?

# Predictor matrices
predmat_csh <- predmat_subdist <- make.predictorMatrix(dat_single_imp)
predmat_csh[, setdiff(
  colnames(predmat_csh), c(sm_predictors, "D1", "D2", "H1", "H2")
)] <- 0
predmat_subdist[, setdiff(
  colnames(predmat_subdist), c(sm_predictors, "newevent", "H1_subdist")
)] <- 0

# A couple of stress tests: (may be a bit silly)
predmat_stress1 <- predmat_stress2 <- make.predictorMatrix(dat_single_imp)

# Cox for rel, ignore NRM (which is our outcome loool)
# Im sick of everything lets ignore outcome!!
predmat_stress1[, setdiff(
  colnames(predmat_stress1), c("PATSEX")#, "D1", "H1")
)] <- 0

predmat_stress1[, setdiff(
  colnames(predmat_stress1), c(sm_predictors)#, "D1", "H1")
)] <- 0

# Cox for NRM (ignoring competing event)
# Now not the hazard, just the time
predmat_stress2[, setdiff(
  colnames(predmat_stress2), c(sm_predictors, "D2", "time_ci_adm")
)] <- 0

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
    return(dato)
  }
)

# First actually check the variance in the subdist marginal haz
rbindlist(cens_imp_dats, idcol = "imp") |>
  ggplot(aes(newtimes, H1_subdist, group = factor(imp), col = factor(imp))) +
  geom_step() +
  theme(legend.position = "none")

# Hardly any difference! Just compute once?

# Now start the imps..
plan(multisession, workers = 3)

# First: imps csh (shall we use the bind mids thing? Combine multiple smcfcs?)
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
    imps$impDatasets[[1]]
  },
  future.seed = TRUE
)

impdats_smcfcs_csh <- future_lapply(
  X = cens_imp_dats,
  FUN = function(imp_times) {
    #imp_times$newevent <- as.numeric(imp_times$newevent) - 1L
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
    imps$impDatasets[[1]]
  },
  future.seed = TRUE
)


# Now the stress tests
impdats_mice_stress1 <- future_lapply(
  X = cens_imp_dats,
  FUN = function(imp_times) {
    imps <- mice(
      data = imp_times,
      m = n_imps,
      maxit = n_cycles,
      method = meths_mice,
      predictorMatrix = predmat_stress1,
      printFlag = FALSE
    )
    complete(imps, "all")[[1]]
  },
  future.seed = TRUE
)

impdats_mice_stress2 <- future_lapply(
  X = cens_imp_dats,
  FUN = function(imp_times) {
    imps <- mice(
      data = imp_times,
      m = n_imps,
      maxit = n_cycles,
      method = meths_mice,
      predictorMatrix = predmat_stress2,
      printFlag = FALSE
    )
    complete(imps, "all")[[1]]
  },
  future.seed = TRUE
)

plan(sequential)



# Coefficients ------------------------------------------------------------


mods_smcfcs_subdist <- lapply(
  impdats_smcfcs_subdist,
  function(impdat) {
    mod <- do.call(
      what = survival::coxph,
      args = list(
        formula = update(sm_form, Surv(newtimes, newevent) ~ .),
        data = impdat,
        x = TRUE
      )
    )
    # Cens times are different in each, so predict at horizon
    #cbind.data.frame(preds = drop(predictRisk(mod, newdata = impdat, times = horizon)), impdat)
  }
)

mods_smcfcs_csh <- lapply(
  impdats_smcfcs_csh,
  function(impdat) {
    mod <- do.call(
      what = survival::coxph,
      args = list(
        formula = update(sm_form, Surv(newtimes, newevent == 2) ~ .),
        data = impdat,
        x = TRUE
      )
    )
    # Cens times are different in each, so predict at horizon
    #cbind.data.frame(preds = drop(predictRisk(mod, newdata = impdat, times = horizon)), impdat)
  }
)

mods_mice_csh <- lapply(
  impdats_mice_csh,
  function(impdat) {
    mod <- do.call(
      what = survival::coxph,
      args = list(
        formula = update(sm_form, Surv(newtimes, newevent == 2) ~ .),
        data = impdat,
        x = TRUE
      )
    )
    # Cens times are different in each, so predict at horizon
    #cbind.data.frame(preds = drop(predictRisk(mod, newdata = impdat, times = horizon)), impdat)
  }
)

mods_mice_subdist <- lapply(
  impdats_mice_subdist,
  function(impdat) {
    mod <- do.call(
      what = survival::coxph,
      args = list(
        formula = update(sm_form, Surv(newtimes, newevent == 2) ~ .),
        data = impdat,
        x = TRUE
      )
    )
    # Cens times are different in each, so predict at horizon
    #cbind.data.frame(preds = drop(predictRisk(mod, newdata = impdat, times = horizon)), impdat)
  }
)

mods_cca <- lapply(
  cens_imp_dats,
  function(impdat) {
    mod <- do.call(
      what = survival::coxph,
      args = list(
        formula = update(sm_form, Surv(newtimes, newevent == 2) ~ .),
        data = impdat,
        x = TRUE
      )
    )
    # Cens times are different in each, so predict at horizon
    #cbind.data.frame(preds = drop(predictRisk(mod, newdata = impdat, times = horizon)), impdat)
  }
)


# Also CCA
rbind(
  cbind(
    broom::tidy(pool(mods_smcfcs_csh), conf.int = TRUE, exponentiate = TRUE),
    "method" = "SMC-FCS CSH"
  ),
  cbind(
    broom::tidy(pool(mods_smcfcs_subdist), conf.int = TRUE, exponentiate = TRUE),
    "method" = "SMC-FCS Subdist"
  ),
  cbind(
    broom::tidy(pool(mods_cca), conf.int = TRUE, exponentiate = TRUE),
    "method" = "CCA"
  ),
  cbind(
    broom::tidy(pool(mods_mice_csh), conf.int = TRUE, exponentiate = TRUE),
    "method" = "MICE CSH"
  ),
  cbind(
    broom::tidy(pool(mods_mice_subdist), conf.int = TRUE, exponentiate = TRUE),
    "method" = "MICE Subdist"
  )
) |>
  ggplot(aes(estimate, term, col = method)) +
  geom_point(size = 3, position = position_dodge(width = 0.75)) +
  geom_linerange(
    linewidth = 1,
    aes(xmin = conf.low, xmax = conf.high),
    position = position_dodge(width = 0.75)
  ) +
  scale_x_continuous(
    name = "Hazard ratio (95% CI)",
    trans = "log",
    breaks = c(0.5, 1, 1.5, 2, 3)
  ) +
  scale_color_brewer(palette = "Dark2") +
  theme_minimal() +
  coord_cartesian(xlim = c(0.65, 3)) +
  geom_vline(xintercept = 1, linetype = "dotted")



# Now the stress test -----------------------------------------------------


mods_mice_stress1 <- lapply(
  impdats_mice_stress1,
  function(impdat) {
    mod <- do.call(
      what = survival::coxph,
      args = list(
        formula = update(sm_form, Surv(newtimes, newevent == 2) ~ .),
        data = impdat,
        x = TRUE
      )
    )
    # Cens times are different in each, so predict at horizon
    #cbind.data.frame(preds = drop(predictRisk(mod, newdata = impdat, times = horizon)), impdat)
  }
)

mods_mice_stress2 <- lapply(
  impdats_mice_stress2,
  function(impdat) {
    mod <- do.call(
      what = survival::coxph,
      args = list(
        formula = update(sm_form, Surv(newtimes, newevent == 2) ~ .),
        data = impdat,
        x = TRUE
      )
    )
    # Cens times are different in each, so predict at horizon
    #cbind.data.frame(preds = drop(predictRisk(mod, newdata = impdat, times = horizon)), impdat)
  }
)


rbind(
  cbind(
    broom::tidy(pool(mods_cca), conf.int = TRUE, exponentiate = TRUE),
    "method" = "CCA"
  ),
  cbind(
    broom::tidy(pool(mods_mice_subdist), conf.int = TRUE, exponentiate = TRUE),
    "method" = "MICE Subdist"
  ),
  cbind(
    broom::tidy(pool(mods_mice_stress1), conf.int = TRUE, exponentiate = TRUE),
    "method" = "MICE stress 1"
  ),
  cbind(
    broom::tidy(pool(mods_mice_stress2), conf.int = TRUE, exponentiate = TRUE),
    "method" = "MICE stress 2"
  )
) |>
  ggplot(aes(estimate, term, col = method)) +
  geom_point(size = 3, position = position_dodge(width = 0.75)) +
  geom_linerange(
    linewidth = 1,
    aes(xmin = conf.low, xmax = conf.high),
    position = position_dodge(width = 0.75)
  ) +
  scale_x_continuous(
    name = "Hazard ratio (95% CI)",
    trans = "log",
    breaks = c(0.5, 1, 1.5, 2, 3)
  ) +
  scale_color_brewer(palette = "Dark2") +
  theme_minimal() +
  coord_cartesian(xlim = c(0.65, 3)) +
  geom_vline(xintercept = 1, linetype = "dotted")


library(rstpm2)

fit <- stpm2(
  Surv(newtimes, newevent == 2) ~ 1
  , data=dat_single_imp, df = 4)
fit <- stpm2(
  update(sm_form, Surv(newtimes, newevent == 2) ~ .)
  , data=dat_single_imp, df = 4)

fit <- stpm2(Surv(time_ci_adm, status_ci_adm == 1) ~ 1, data = dat_single_imp, df = 4)
fit <- stpm2(Surv(time_ci_adm, status_ci_adm == 2) ~ 1, data = dat_single_imp, df = 4)
grid <- seq(0.01, 60, by = 0.1)

plot(
  grid,
  predict(
  fit,
  type = "hazard",
  newdata = data.frame(newtimes = grid))
)


plot(
  grid,
  predict(
    fit,
    type = "hazard",
    newdata = data.frame(time_ci_adm = grid))
)


# Predictions -------------------------------------------------------------



# This won't work because of different censoring times
# (when objects were mice)
#Reduce(mice::ibind, imps_mice_csh)


# Try predict tings first with one of them
horizon <- 24 # 2 years

cloglog <- function(x) log(-log(1 - x))
inv_cloglog <- function(x) 1 - exp(-exp(x))

preds_mice_csh <- lapply(
  impdats_mice_csh,
  function(impdat) {
    mod <- do.call(
      what = survival::coxph,
      args = list(
        formula = update(sm_form, Surv(newtimes, newevent == 2) ~ .),
        data = impdat,
        x = TRUE
      )
    )
    # Cens times are different in each, so predict at horizon
    cbind.data.frame(preds = drop(predictRisk(mod, newdata = impdat, times = horizon)), impdat)
  }
)

preds_mice_subdist <- lapply(
  impdats_mice_subdist,
  function(impdat) {
    mod <- do.call(
      what = survival::coxph,
      args = list(
        formula = update(sm_form, Surv(newtimes, newevent == 2) ~ .),
        data = impdat,
        x = TRUE
      )
    )
    # Cens times are different in each, so predict at horizon
    cbind.data.frame(preds = drop(predictRisk(mod, newdata = impdat, times = horizon)), impdat)
  }
)

preds_smcfcs_subdist <- lapply(
  impdats_smcfcs_subdist,
  function(impdat) {
    mod <- do.call(
      what = survival::coxph,
      args = list(
        formula = update(sm_form, Surv(newtimes, newevent) ~ .),
        data = impdat,
        x = TRUE
      )
    )
    # Cens times are different in each, so predict at horizon
    cbind.data.frame(preds = drop(predictRisk(mod, newdata = impdat, times = horizon)), impdat)
  }
)


df_mice_csh <- rbindlist(preds_mice_csh, idcol = "imp")
df_mice_subdist <- rbindlist(preds_mice_subdist, idcol = "imp")
df_smcfcs_subdist <- rbindlist(preds_smcfcs_subdist, idcol = "imp")

test <- cbind.data.frame(
  preds_csh = df_mice_csh[, .(
    preds_pooled = inv_cloglog(mean(cloglog(preds)))
  ), by = id][["preds_pooled"]],
  preds_subdist = df_mice_subdist[, .(
    preds_pooled = inv_cloglog(mean(cloglog(preds)))
  ), by = id][["preds_pooled"]],
  preds_smcfcs_subdist = df_smcfcs_subdist[, .(
    preds_pooled = inv_cloglog(mean(cloglog(preds)))
  ), by = id][["preds_pooled"]]
)

plot(test$preds_smcfcs_subdist, test$preds_subdist)
abline(a = 0, b = 1, col = "blue")

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





# Old ---------------------------------------------------------------------





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
