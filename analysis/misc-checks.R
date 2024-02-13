# All this is separate qmd
source(here::here("packages.R"))
library(lme4)
library(GLMMadaptive)
library(glmmTMB)
library(rms)
invisible(lapply(list.files(here::here("R"), full.names = TRUE), source))
tar_load(applied_dat_raw)

# Overall missing vals
naniar::gg_miss_var(applied_dat_raw, show_pct = TRUE)
print(naniar::miss_var_summary(applied_dat_raw), n = 60)


# CCA ind -----------------------------------------------------------------

#
#https://stats.stackexchange.com/questions/398869/r-how-to-fit-a-glmm-in-nlme
# try also brm
preds_orig <- c(
  "hctci_risk",
  "age_allo1_decades", # Already centered at age = 60 (median age 58)
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

preds_extra <- c(
  "agedonor_allo1_decades",
  "vchromos_preallo1"
  "DONSEX_allo1_1"
  # Don't show the next two in the summary:
  "year_allo1",
  "intdiagtr_allo1"
  #submps_allo1, tbi_allo1, source_allo1 not included, whatever
)

tar_load(applied_dat_raw)

applied_dat_raw[, tceldepl_bin := ifelse(tceldepl_allo1 == "no", "no", "yes")]

par(mfrow = c(1, 2))
plot(
  prodlim(
    Hist(time_ci_adm, status_ci_adm) ~ hctci_risk,
    data = applied_dat_raw
  ),
  atrisk = FALSE,
  legend.cex = 0.7,
  cause = 1
)

coxph(
  Surv(time_ci_adm, status_ci_adm == 1) ~ age_allo1_decades +
    PATSEX +
    hctci_risk +
    bmi,
  data = applied_dat_raw
)

library(rms)
dd <- datadist(applied_dat_raw)
options(datadist = "dd")
dd$limits$bmi[2] <- 22

mod <- cph(
  Surv(time_ci_adm, status_ci_adm == 1) ~ age_allo1_decades +
    PATSEX +
    #hctci_risk +
    #WEIGLOSS_allo1 +
    rcs(bmi, 4),
  data = applied_dat_raw
)

# Bmi had effect on relapse>

mod
ggplot(Predict(mod, bmi, ref.zero = TRUE))

# They got more gvhd.. and they are more overweight?, and relatively
# more high DIPSS?
# many  more splenecotimized in hctci_high?
par(mfrow = c(1, 2))
plot(
  prodlim(
    Hist(time_ci_agvh2to4, status_ci_agvh2to4) ~ hctci_risk,
    data = applied_dat_raw
  ),
  atrisk = FALSE,
  legend.cex = 0.7,
  cause = 1
)

cmprsk::cuminc(applied_dat_raw$time_ci_agvh2to4, applied_dat_raw$status_ci_agvh2to4,
               applied_dat_raw$hctci_risk)$Tests |> round(digits = 4)

library(gtsummary)
applied_dat_raw[, -c("id_new", "CENTRE_allo1", "AA_CTY")] |>
  tbl_summary(by = "hctci_risk") |>
  add_p()

cmprsk::cuminc(applied_dat_raw$time_ci_adm, applied_dat_raw$status_ci_adm,
               applied_dat_raw$hctci_risk)$Tests |> round(digits = 4)

library(gtsummary)
applied_dat_raw[, c("sweat_allo1", "const_allo1", "WEIGLOSS_allo1")] |> #, "pb_allo1")] |>
  tbl_summary(by = "const_allo1") #|>
 # add_p()

tar_load(applied_dat)
sm_predictors <- applied_dat$sm_predictors
df <- applied_dat$dat
applied_dat_raw[, -c("id_new", "CENTRE_allo1", "AA_CTY")] |>
  tbl_summary(by = "DONSEX_allo1_1") |>
  add_p()

coxph(
  Surv(time_ci_adm, status_ci_adm == 2) ~ I(plat_allo1  > 500) +
    pb_allo1 +
    age_allo1_decades,
  data = applied_dat_raw
)

library(splines)
coxph(
  Surv(time_ci_adm, status_ci_adm == 1) ~
    sweat_allo1 *
    ns(pb_allo1 + 0.1, 3),
  data = applied_dat_raw
)

library(rms)
d
mod <- cph(
  Surv(time_ci_adm, status_ci_adm == 1) ~
    sweat_allo1 +
    pb_allo1,
    #rcs(pb_allo1, 4),
  data = applied_dat_raw
)




applied_dat_raw[!is.na(DONSEX_allo1_1)] |>
  ggplot(aes(agedonor_allo1_decades, col = DONSEX_allo1_1)) +
  geom_density() #+
  #facet_wrap(~ sweat_allo1)

applied_dat_raw[!is.na(DONSEX_allo1_1)] |>
  ggplot(aes(agedonor_allo1_decades, col = donrel_bin)) +
  geom_density()


ggplot(Predict(mod, pb_allo1, const_allo1, ref.zero = TRUE))

preds <- c(
  "age_allo1_decades",
  "agedonor_allo1_decades",
  "year_allo1",
  "intdiagtr_allo1",
  "tbi_allo1",
  "ric_allo1",
  "donrel_bin",
  "PATSEX",
  "DONSEX_allo1_1",
  "vchromos_preallo1",
  "WEIGLOSS_allo1",
  "hb_allo1",
  "pb_allo1",
  "wbc_allo1",
  "KARNOFSK_threecat",
  "cmv_match",
  "submps_allo1"
)
applied_dat_raw$CCA_ind <- as.numeric(complete.cases(applied_dat_raw[, ..preds]))
mean(applied_dat_raw$CCA_ind)

dd <- datadist(applied_dat_raw)
options(datadist = "dd")
dd$limits$intdiagtr_allo1[2] <- min(applied_dat_raw$intdiagtr_allo1)

mod <- lrm(
  CCA_ind ~ rcs(year_allo1, 4) +
    log(intdiagtr_allo1 + 0.01) +
    status_os_adm *
    rcs(time_os_adm, 4) +
    submps_allo1 +
    donrel_bin +
    KARNOFSK_threecat +
    ric_allo1,
  data = applied_dat_raw
)

mod

theme_set(theme_light(base_size = 14))
ggplot(Predict(mod, intdiagtr_allo1, ref.zero = TRUE), ylim = c(-2.5, 2.5)) +
  geom_hline(yintercept = 0, linetype = "dotted")

ggplot(Predict(mod, year_allo1, ref.zero = TRUE), ylim = c(-2.5, 2.5)) +
  geom_hline(yintercept = 0, linetype = "dotted")

ggplot(Predict(mod, time_os_adm, status_os_adm, ref.zero = TRUE), ylim = c(-2.5, 2.5)) +
  geom_hline(yintercept = 0, linetype = "dotted")

anova(mod)


# Importance of donor age
mod_donrel <- cph(
  Surv(time_ci_adm, status_ci_adm == 1) ~ donrel_bin +
    DONSEX_allo1_1 +
    #agedonor_allo1_decades +
    age_allo1_decades +
    PATSEX +
    cmv_match +
    ric_allo1 +
    KARNOFSK_threecat,
  data = applied_dat_raw
)
mod_donrel


mod_flex_donrel <- cph(
  Surv(time_ci_adm, status_ci_adm == 1) ~ donrel_bin *
    DONSEX_allo1_1 *
    agedonor_allo1_decades +
    #rcs(agedonor_allo1_decades, 4) +
    age_allo1_decades +
    PATSEX, # +
    #cmv_match +
    #ric_allo1 +
    #KARNOFSK_threecat,
  data = applied_dat_raw
)
mod_flex_donrel
ggplot(
  Predict(
    mod_flex_donrel,
    agedonor_allo1_decades,
    DONSEX_allo1_1,
    donrel_bin,
    ref.zero = TRUE
  )
)

mod_miss <- glmmTMB(
  !is.na(hctci_risk) ~ ns(year_allo1, 4) +
    PATSEX +
    intdiagtr_allo1 +
    donrel_bin +
    ric_allo1 +
    cmv_match +
    submps_allo1 +
    tbi_allo1 +
    KARNOFSK_threecat +
    status_os_adm * time_os_adm +
    (1 | AA_CTY / CENTRE_allo1), #/ CENTRE_allo1),
  data = applied_dat_raw,
  family = binomial(link = "logit")
)
summary(mod_miss)
performance::r2_nakagawa(mod_miss)

newdat <- expand.grid(
  year_allo1 = 10,
  status_os_adm = c(0, 1),
  AA_CTY = unique(applied_dat_raw$AA_CTY),
  time_os_adm = 0:60
)

cbind.data.frame(
  newdat,
  preds = predict(mod_miss, newdata = newdat),
  preds_marg = predict(mod_miss, newdata = newdat, re.form = NA)
) |>
  ggplot(aes(time_os_adm, preds)) +
  geom_line(
    aes(
      group = interaction(AA_CTY, status_os_adm),
      linetype = factor(status_os_adm)
    ),
    linewidth = 1, alpha = 0.25
  ) +
  geom_line(
    aes(
      y = preds_marg,
      linetype = factor(status_os_adm),
      col =  factor(status_os_adm)
    ),
    linewidth = 3
  ) +
  facet_wrap(~ status_os_adm)



# Now the glmer
mod_miss <- glmer(
  #CCA_ind ~
  as.numeric(is.na(hctci_risk)) ~
    ns(year_allo1, 4) +
    intdiagtr_allo1 +
    #tbi_allo1 +
    factor(status_os_adm) +
    time_os_adm +
    #(1 | AA_CTY),
    + (1 | AA_CTY / CENTRE_allo1),
  data = applied_dat_raw,
  family = binomial
)
summary(mod_miss)
lattice::dotplot(ranef(mod_miss))
performance::r2_nakagawa(mod_miss)

library(splines)

mod_miss <- glmer(
  CCA_ind ~
  #as.numeric(is.na(hctci_risk)) ~
    year_allo1 +
    factor(status_os_adm) *
    ns(time_os_adm, 4) + (ns(time_os_adm, 4) | AA_CTY),
    #(time_os_adm | AA_CTY),
   # + (year_allo1 | AA_CTY / CENTRE_allo1), #+
  data = applied_dat_raw,
  family = binomial
)
summary(mod_miss)
newdat <- expand.grid(
  year_allo1 = 10,
  status_os_adm = c(0, 1),
  AA_CTY = unique(applied_dat_raw$AA_CTY),
  time_os_adm = 0:60
)

cbind.data.frame(
  newdat,
  preds = predict(mod_miss, newdata = newdat),
  preds_marg = predict(mod_miss, newdata = newdat, re.form = NA)
) |>
  ggplot(aes(time_os_adm, preds)) +
  geom_line(
    aes(group = interaction(AA_CTY, status_os_adm),
        linetype = factor(status_os_adm)),
    linewidth = 1, alpha = 0.25
  ) +
  geom_line(aes(y = preds_marg, linetype = factor(status_os_adm),
                col =  factor(status_os_adm)),
            linewidth = 3)


lattice::dotplot(ranef(mod_miss))
performance::r2_nakagawa(mod_miss)



# Other (to sort) ---------------------------------------------------------



coxph(
  Surv(time_ci_adm, status_ci_adm == 1) ~ relevel(
    interaction(
      DONSEX_allo1_1, PATSEX
    ), ref = "Male.Male"
  ) +
    age_allo1_decades,
  data = applied_dat_raw[donrel_bin == "Other"]
)

# Add: cuminc curves per var?
# Proportionality check?
library(lme4)
library(GLMMadaptive)



library(rms)

dd <- datadist(applied_dat_raw)
options(datadist="dd")
dd$limits$intdiagtr_allo1[2] <- min(applied_dat_raw$intdiagtr_allo1)

mod <- lrm(
  CCA_ind ~ rcs(year_allo1, 4) +
    log(intdiagtr_allo1 + 0.01) +
    status_os_adm *
    rcs(time_os_adm, 4) +
    rcs(year_allo1, 4) +
    submps_allo1 +
    donrel_bin +
    KARNOFSK_threecat +
    ric_allo1,
  data = applied_dat_raw
)

mod

theme_set(theme_light(base_size = 14))
ggplot(Predict(mod, intdiagtr_allo1, ref.zero = TRUE), ylim = c(-2.5, 2.5)) +
  geom_hline(yintercept = 0, linetype = "dotted")
ggplot(Predict(mod, year_allo1, ref.zero = TRUE), ylim = c(-2.5, 2.5)) +
  geom_hline(yintercept = 0, linetype = "dotted")
ggplot(Predict(mod, time_os_adm, status_os_adm, ref.zero = TRUE), ylim = c(-2.5, 2.5)) +
  geom_hline(yintercept = 0, linetype = "dotted")
anova(mod)

mod_miss <- mixed_model(
  fixed = #as.numeric(is.na(hb_allo1)) ~ year_allo1 +
    CCA_ind ~ year_allo1 +
    intdiagtr_allo1 +
    age_allo1_decades +
    tbi_allo1 +
    ric_allo1 +
    donrel_bin +
    PATSEX +
    cmv_match +
    relevel(factor(status_os_adm), ref = "1") *
    time_os_adm,
  random = ~ 1 | AA_CTY,
  data = applied_dat_raw,
  family = binomial()
)
summary(mod_miss)
lattice::dotplot(ranef(mod_miss)[order(ranef(mod_miss)), ])

newdat <- data.frame("year_allo1" = 0, "CENTRE_allo1" = unique(applied_dat_raw$CENTRE_allo1))
plot(newdat$CENTRE_allo1, predict(mod_miss, ))

summary(mod_miss)


mod_miss <- glmer(
  as.numeric(is.na(hctci_risk)) ~ year_allo1 +
    #intdiagtr_allo1 +
    #tbi_allo1 +
    #factor(status_os_adm) +
    #time_os_adm +
     #(1 | AA_CTY),
    + (1 | AA_CTY / CENTRE_allo1),
  data = applied_dat_raw,
  family = binomial
)
lattice::dotplot(ranef(mod_miss))
summary(mod_miss)

performance::r2_nakagawa(mod_miss)

partR2::partR2(
  mod_miss,
  R2_type = "conditional",
  partvars = c("intdiagtr_allo1", "time_os_adm")
)

class(mod_miss)
summary(mod_miss)
lattice::dotplot(ranef(mod_miss))


# glmmtmb?
mod_miss <- glmer(
  as.numeric(is.na(hb_allo1)) ~ year_allo1 +
    intdiagtr_allo1 +
    age_allo1_decades +
    #donrel_bin +
    PATSEX +
    #cmv_match +
    #tbi_allo1 +
    #ric_allo1 +
    factor(status_ci_adm) +
    time_ci_adm +
    (1 | CENTRE_allo1),
  data = applied_dat_raw,
  family = binomial
)
# Do it by country??
summary(mod_miss)
names(applied_dat_raw)

lattice::dotplot(ranef(mod_miss))
library(sjPlot)
plot_model(mod_miss, type = "re", transform = "plogis")


summary(mod_miss)


# Test frailty ------------------------------------------------------------


library(JointAI)

dat_test <- applied_dat_raw[, c(
  "time_os_adm",
  "status_os_adm",
  "AA_CTY",
  "age_allo1_decades",
  "PATSEX",
  "intdiagtr_allo1",
  "year_allo1",
  "agedonor_allo1_decades"
)] |> data.frame()

future::plan(future::multisession, workers = 3)
test1 <- coxph_imp(
  Surv(time_os_adm, status_os_adm) ~ agedonor_allo1_decades +
    age_allo1_decades +
    PATSEX +
    intdiagtr_allo1 +
    year_allo1,
  data = dat_test,
  n.iter = 100,
  n.chains = 3
)
future::plan(future::sequential)

summary(test1)
traceplot(test1)
densplot(test1)


future::plan(future::multisession, workers = 3)
test2 <- coxph_imp(
  Surv(time_os_adm, status_os_adm) ~ agedonor_allo1_decades +
    age_allo1_decades +
    PATSEX +
    intdiagtr_allo1 +
    year_allo1 +
    (1 | AA_CTY),
  data = dat_test,
  n.adapt = 500,
  n.iter = 500,
  n.chains = 3
)
future::plan(future::sequential)

traceplot(test2)
densplot(test1)
summary(test2)
