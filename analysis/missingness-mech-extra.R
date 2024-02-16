# Random effects models + other small ideas to keep for later
source(here::here("packages.R"))
library(lme4)
library(GLMMadaptive)
library(glmmTMB)
library(rms)
library(brms)
invisible(lapply(list.files(here::here("R"), full.names = TRUE), source))
options(contrasts = rep("contr.treatment", 2))

cols <- c("#CABEE9", "#7C7189", "#FAE093", "#D04E59", "#BC8E7D", "#2F3D70")
theme_set( # set base size to 14 instead of 16
  theme_light(base_size = 16, base_family = "Roboto Condensed") +
    theme(
      strip.background = element_rect(fill = cols[2], colour = "white"),
      strip.text = element_text(colour = 'white')
    )
)

tar_load(applied_dat)
tar_load(applied_dat_raw)
dat <- applied_dat$dat



# The donor age ting ------------------------------------------------------

applied_dat_raw[!is.na(donrel_bin) & !is.na(DONSEX_allo1_1)] |>
  ggplot(aes(agedonor_allo1_decades, fill = donrel_bin)) +
  geom_density(alpha = 0.5)

applied_dat_raw[!is.na(donrel_bin) & !is.na(DONSEX_allo1_1)] |>
  ggplot(aes(agedonor_allo1_decades, fill = DONSEX_allo1_1)) +
  geom_density(alpha = 0.5)


# Exploratory model for censoring? (about times imp) ----------------------

mod_cens <- coxph(
  Surv(time_ci_adm, status_ci_adm == 0) ~ year_allo1_decades +
    intdiagallo_decades +
    donrel_bin +
    age_allo1_decades +
    PATSEX +
    cmv_match +
    ric_allo1 +
    KARNOFSK_threecat,
  data = applied_dat$dat
)
summary(mod_cens)

applied_dat$dat[, year_allo1_fac := factor(
  year_allo1_decades * 10,
  labels = 2009:2019
)]

sf_cens <- prodlim(
  Hist(time_ci_adm, status_ci_adm == 0) ~ year_allo1_fac,
  data = applied_dat$dat
)
plot(
  sf_cens,
  limit = 15,
  atrisk = FALSE,
  lty = 1:10,
  legend.x = "bottomleft",
  legend.cex = 0.65
)



# Interesting things ------------------------------------------------------


#https://stats.stackexchange.com/questions/398869/r-how-to-fit-a-glmm-in-nlme

preds <- applied_dat$sm_predictors
dat$CCA_ind <- as.numeric(complete.cases(dat[, ..preds]))


mod_miss <- brm(
  as.numeric(!is.na(hctci_risk)) ~ #ns(year_allo1, 4) +
    #PATSEX +
    #intdiagtr_allo1 +
    #donrel_bin +
    #ric_allo1 +
    #cmv_match +
    #submps_allo1 +
    #tbi_allo1 +
    #KARNOFSK_threecat +
    #status_os_adm * time_os_adm +
    year_allo1 +
    (1 | AA_CTY / CENTRE_allo1), #/ CENTRE_allo1),
  data = data.frame(applied_dat_raw),
  family = bernoulli
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
