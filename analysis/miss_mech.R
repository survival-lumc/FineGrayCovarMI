dat <- generate_dataset(
  n = 2500,
  args_event_times = list(
    mechanism = "correct_FG",
    params = params,
    censoring_type = "none"#,
    #censoring_params = list("exponential" = "0.1 * exp(1 * (as.numeric(X) - 1))")
  ),
  args_missingness = list(
    mech_params = list("prob_missing" = 0.4, "mechanism_expr" = "1.5 * Z")
  )
)

hist(dat$time, breaks = 100, xlim = c(0, 10))
mu <- mean(log(dat$time))
sig <- sd(log(dat$time))
exp(mu + seq(-4, 2) * sig)

hist(scale(log(dat$time)), breaks = 100)
hist(scale(log(dat$time)), breaks = 100)


hist(scale(log10(dat$time)), breaks = 100)

hist(scale(dat$time))
hist(sqrt(dat$time))






# Test Polverelli ---------------------------------------------------------


dat_raw <- tar_read(applied_dat_raw)[complete.cases(status_ci_adm)]
applied_tings <- tar_read(applied_dat)
dat <- process_applied_dat(dat_raw)$dat
dat[, "centre" := dat_raw$CENTRE_allo1]
dat[, "centre_lab" := as.character(as.numeric(centre))]

sm_predictors <- applied_tings$sm_predictors


dat$CCA_ind <- as.numeric(complete.cases(dat))
dat$status_ci_adm_bis <- relevel(factor(dat$status_ci_adm), ref = "2")
prop.table(table(dat$status_ci_adm))
naniar::gg_miss_var(dat, show_pct = TRUE)

# Model for the censoring
library(ggeffects)
library(splines)
library(glmmTMB)

options(contrasts = rep("contr.treatment", 2))

aj_miss <- prodlim(Hist(time_ci_adm, status_ci_adm == 0) ~ CCA_ind + factor(year_allo1_decades), data = dat)
plot(aj_miss, legend.cex = 0.75)
#plot(aj_miss, cause = 1, legend.cex = 0.75)

yearz <- sort(unique(dat$year_allo1_decades))
for (i in seq_along(yearz)) {
  aj_miss <- prodlim(
    Hist(time_ci_adm, status_ci_adm) ~ CCA_ind,
    data = dat[year_allo1_decades == yearz[i]]
  )
  #plot(aj_miss, legend.cex = 0.75)
  plot(aj_miss, cause = 2, legend.cex = 0.75)
  title(yearz[i] * 10 + 2019)
}

mod_cca <- glm(
  CCA_ind ~ year_allo1_decades * status_ci_adm_bis * time_ci_adm + #ns(time_ci_adm, 3) + #* ns(time_ci_adm, 3) +
    intdiagallo_decades + donrel_bin + submps_allo1,# +
    #PATSEX + age_allo1_decades +
    #ric_allo1 + tceldepl_bin + cmv_match + KARNOFSK_threecat,
  data = dat,
  family = binomial()
)
summary(mod_cca)
ggeffects::predict_response(
  mod_cca,
  terms = c("time_ci_adm",
            "status_ci_adm_bis", "year_allo1_decades[-1, -0.5, 0]"),
  margin = "mean_mode"
) |>
  plot()

# Try now with ranef
mod_CCA_centre <-  glmmTMB(
  CCA_ind ~ year_allo1_decades + status_ci_adm * ns(time_ci_adm, 3) +
    intdiagallo_decades + donrel_bin + submps_allo1 +
    PATSEX + age_allo1_decades +
    ric_allo1 + tceldepl_bin + cmv_match + KARNOFSK_threecat +
    (1 | centre), # can also do (1 | country / center) for country effects, with vcenland var
  data = dat,
  family = binomial()
)
summary(mod_CCA_centre)


mod_cens <-




mod_CCA_centre <-  glmmTMB(
  CCA_ind ~ year_allo1_decades + status_ci_adm * ns(time_ci_adm, 4) +
    donrel_bin + tceldepl_binyes
    (1 | centre_lab), # can also do (1 | country / center) for country effects, with vcenland var
  data = dat,
  family = binomial()
)
lme4:::dotplot.ranef.mer(ranef(mod_CCA_centre)$cond)

test <- ggeffects::predict_response(
  mod_CCA_centre,
  terms = c("time_ci_adm [all]", "status_ci_adm"),
  #type = "simulate"
 # margin = "empirical"#,
  type = "fixed"
)

plot(test)

summary(mod_CCA_centre)
summary(mod_cca)

mod_CCA_centre <-  glmmTMB(
  CCA_ind ~ year_allo1_decades +
    intdiagallo_decades +
    (1 | centre_lab), # can also do (1 | country / center) for country effects, with vcenland var
  data = dat,
  family = binomial()
)
plot(residuals(mod_CCA_centre))

library(DHARMa)
res1 <- simulateResiduals(mod_CCA_centre)
plot(res1)
lme4:::dotplot.ranef.mer(ranef(mod_CCA_centre)$cond)

#
ranef(mod_CCA_centre)$cond

dat[, .(
  size_centre = .N,
  prop_complete = mean(complete.cases(.SD))
), by = centre_lab][order(size_centre)]

dat[, .(
  size_centre = .N,
  prop_complete = mean(complete.cases(.SD))
), by = centre] |>
  ggplot(aes(qlogis(prop_complete), reorder(centre, prop_complete))) +
  geom_point(aes(size = size_centre))


library(ggrepel)
dat[, .(
  size_centre = .N,
  prop_complete = mean(complete.cases(.SD))
), by = centre] |>
  ggplot(aes(size_centre, reorder(centre, size_centre))) +
  geom_point(aes(size = prop_complete), alpha = 0.25) +
  geom_text_repel(
    aes(label = centre),
    size = 2,
    max.overlaps = 25
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank()
  ) +
  labs(y = NULL)

dat[, .(
  size_centre = .N,
  prop_complete = mean(complete.cases(.SD))
), by = centre] |>
  ggplot(aes(size_centre, prop_complete)) +
  geom_point(aes(size = size_centre), alpha = 0.25) +
  geom_text_repel(
    aes(label = centre),
    size = 2,
    max.overlaps = 25
  ) +
  coord_cartesian(ylim = c(-0.5, 1.5), xlim = c(-50, 225))
  theme_minimal() +
  theme(
    #axis.text.y = element_blank()
  ) #+
  #labs(y = NULL)


dat[, .(
    size_centre = .N,
    prop_complete = mean(complete.cases(.SD))
  ), by = centre] |>
  ggplot(aes(log(size_centre), prop_complete)) +
  geom_point() +
  geom_smooth()


ggeffects::predict_response(
  mod_CCA_centre,
  terms = c("time_ci_adm [all]", "status_ci_adm"), margin = "mean_mode"
) |>
  plot()




form_cens <- reformulate(
  response = "Surv(time_ci_adm, status_ci_adm == 0)",
  termlabels = sm_predictors
)


# Modle for the cens:
mod_cens <- coxph(
  form_cens,
  data = dat
)
summary(mod_cens)



dat_splitmonth <- survSplit(
  Surv(time_ci_adm, status_ci_adm == 0) ~ CCA_ind +
    year_allo1_decades +
    intdiagallo_decades + donrel_bin + submps_allo1 +
    PATSEX + age_allo1_decades +
    ric_allo1 + tceldepl_bin + cmv_match + KARNOFSK_threecat + centre_lab,
  data = dat,
  cut = seq(1, 59, by = 1),
  episode = "month"
)

modo <- coxph(
  Surv(tstart, time_ci_adm, event) ~ CCA_ind + CCA_ind:month +
    year_allo1_decades +
    intdiagallo_decades + donrel_bin + submps_allo1 +
    PATSEX + age_allo1_decades +
    ric_allo1 + tceldepl_bin + cmv_match + KARNOFSK_threecat +
    frailty.gaussian(centre_lab),
  data = dat_splitmonth
)
summary(modo)
predict_response(modo, terms = c("CCA_ind", "month")) |>
  plot()

dat_test <- copy(dat)
dat_test$stat_ci <- as.numeric(as.character(dat_test$status_ci_adm))

hi <- prodlim(Hist(time_ci_adm, stat_ci) ~ donrel_bin, data = dat_test) #|>
plot(hi, cause = 1)
plot(hi, cause = 2)

coxph(Surv(time_ci_adm, stat_ci == 2) ~ donrel_bin, data = dat_test)
coxph(Surv(time_ci_adm, stat_ci == 1) ~ donrel_bin, data = dat_test)

# Now for cens
dat_splitmonth <- survSplit(
  Surv(time_ci_adm, status_ci_adm == 2) ~ CCA_ind +
    year_allo1_decades +
    intdiagallo_decades + donrel_bin + submps_allo1 +
    PATSEX + age_allo1_decades +
    ric_allo1 + tceldepl_bin + cmv_match + KARNOFSK_threecat,
  data = dat,
  cut = seq(1, 59, by = 1),
  episode = "month"
)

modo <- coxph(
  Surv(tstart, time_ci_adm, event) ~ CCA_ind  +
    CCA_ind:nsk(month, df = 4) + year_allo1_decades, # +
    #year_allo1_decades +
    #intdiagallo_decades + donrel_bin + submps_allo1 +
    #PATSEX + age_allo1_decades +
    #ric_allo1 + tceldepl_bin + cmv_match + KARNOFSK_threecat,
  data = dat_splitmonth
)
summary(modo)

tdata <- expand.grid(
  CCA_ind = 1,
  times = seq(0, 60, by = 1),
  year_allo1_decades = 0
)

tdata$month <- tdata$times

yhat <- predict(modo, newdata = tdata, se.fit = TRUE, reference = "zero")
yy <- yhat$fit + outer(yhat$se.fit, c(0, -1.96, 1.96), '*')

df_hrs_step <- data.table(yy, "times" = tdata$times, CCA_ind = tdata$CCA_ind)
setnames(df_hrs_step, c("log_hr", "lower", "upper", "times", "CCA_ind"))


df_hrs_step |>
  ggplot(aes(times, log_hr, group = CCA_ind)) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    fill = "blue",
    alpha = 0.25,
    col = NA
  ) +
  geom_line(linewidth = 1.25, col = "blue") +
  scale_y_continuous(
    breaks = log(c(0.5, 0.75, 1, 1.25, 1.5, 2, 3)),
    labels = c(0.5, 0.75, 1, 1.25, 1.5, 2, 3)
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "Time since alloSCT (months)",
    y = "Hazard ratio (95% CI) Class II vs 10 out 10"
  ) +
  #coord_cartesian(ylim = log(c(0.4, 3))) +
  ylab(NULL) +
  xlab(NULL)
