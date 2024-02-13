# Make into quarto

# Checking how many imputations
source(here::here("packages.R"))
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

impdats <- data.table(tar_read(applied_impdats))
tar_load(applied_dat)
dat <- applied_dat$dat
impdats[, .(.N), by = c("tar_batch", "tar_rep", "method")]

#tar_meta() |> View() # checking warnings


# Non-param curves --------------------------------------------------------

np_curves <- prodlim(
  Hist(time_ci_adm, status_ci_adm) ~ 1,
  data = dat
)

plot(np_curves, cause = "stacked")
plot(np_curves, cause = 1, col = "blue") # note the jump in cont prog pats
plot(np_curves, cause = 2, add = TRUE, lty = 2)


# Descriptives table? -----------------------------------------------------


naniar::gg_miss_upset(dat)
naniar::miss_var_summary(dat)
naniar::prop_complete_case(dat)

sm_predictors <- applied_dat$sm_predictors
df_predictors <- data.frame(dat[, ..sm_predictors])
miss_props <- sapply(df_predictors, function(col) mean(is.na(col)))
df_predictors[, names(miss_props)[order(miss_props, decreasing = TRUE)]] |>
  naniar::vis_miss()

ggmice::plot_pattern(df_predictors, npat = 100, rotate = TRUE)



# Check imputed values vs observed? For continuous? -----------------------




# Pooled tings ------------------------------------------------------------




sm_form_fg <- reformulate(
  termlabels = applied_dat$sm_predictors,
  response = "Surv(newtimes, newevent)"
)
sm_form_cs1 <- update(sm_form_fg, Surv(time_ci_adm, status_ci_adm == 1) ~ .)
sm_form_cs2 <- update(sm_form_fg, Surv(time_ci_adm, status_ci_adm == 2) ~ .)

mods_imp_dats <- impdats[, .(
  mods_fg = list(coxph(sm_form_fg, data = .SD, x = TRUE)),
  mods_cs1 = list(coxph(sm_form_cs1, data = .SD, x = TRUE)),
  mods_cs2 = list(coxph(sm_form_cs2, data = .SD, x = TRUE))
), by = c("tar_batch", "tar_rep", "method")] |>
  melt(
    measure.vars = c("mods_fg", "mods_cs1", "mods_cs2"),
    variable.name = "mod_type",
    value.name = "mod"
  )

times_preds <- seq(0.01, 60, by = 0.1)
basecuminc_imp_dats <- mods_imp_dats[mod_type == "mods_fg", .(
  "base" = list(
    cbind.data.frame(
      time = times_preds,
      cuminc = 1 - predictCox(object = mod[[1]], times = times_preds, centered = FALSE)$survival
    )
  )
), by = c("tar_batch", "tar_rep", "method")]

basecuminc_imp_dats[, unlist(base, recursive = FALSE), by = c("tar_batch", "tar_rep", "method")][, .(
  pred = inv_cloglog(mean(cloglog(cuminc)))
), by = c("method", "time")] |>
  ggplot(aes(time, pred, group = method, col = method)) +
  geom_step(aes(linetype = method), linewidth = 1) +
  scale_color_manual(values = cols[c(1, 2, 6, 4, 5)]) +
  coord_cartesian(ylim = c(0, 0.25)) +
  labs(
    col = "Method",
    linetype = "Method",
    x = "Time since alloHCT (months)",
    y = "Baseline cumulative incidence function"
  )

# Exponentiate here if needed
pooled_mods <- mods_imp_dats[, .(
  summ = list(tidy(pool(mod), conf.int = TRUE, exponentiate = TRUE))
), by = c("method", "mod_type")][, unlist(summ, recursive = FALSE), by = c("method", "mod_type")]

pooled_mods[!(term %in% c("intdiagallo_decades", "year_allo1_decades")) & mod_type == "mods_fg"] |>
  ggplot(aes(term, estimate, group = method, col = method, shape = method)) +
  geom_point(position = position_dodge(width = 0.75)) +
  geom_pointrange(
    aes(ymin = conf.low, ymax = conf.high),
    position = position_dodge(width = 0.75)
  ) +
  coord_flip() +
  geom_hline(yintercept = 1, linetype = "dotted") +
  scale_color_manual(values = cols[c(1, 2, 6, 4, 5)]) +
  theme(legend.position = "top") +
  scale_x_discrete(limits = rev) +
  scale_y_continuous(
   trans = "log",
   breaks = c(0.5, 1, 1.5, 2, 3)
   #breaks = log(c(0.5, 0.75, 1, 1.5, 2, 3)),
   #labels = c(0.5, 0.75, 1, 1.5, 2, 3)
  )


mods_imp_dats[mod_type == "mods_fg" & method == "SMC-FCS Fine-Gray"]$mod |>
  pool() |>
  howManyImputations::how_many_imputations()


coxph(
  Surv(time_ci_adm, status_ci_adm == 1) ~ hctci_risk   +
    age_allo1_decades +
    donrel_bin +
    PATSEX +
    cmv_match +
    submps_allo1 +
    ric_allo1 +
    KARNOFSK_threecat +
    tceldepl_bin,
  data = applied_dat$dat
)

gtsummary::

library(gtsummary)

applied_dat$dat[, !c("id")] |>
  tbl_summary()

applied_dat$dat[, !c("id")] |>
  tbl_summary(by = "hctci_risk") |>
  add_p()

# Eventually do not show submps,
# in supplement, table with effects + CIs on cause-specific hazards

# Proportionality checks --------------------------------------------------


cox.zph(mods_imp_dats[1, ]$mod[[1]], terms = TRUE)

par(mfrow = c(4, 4))
prop <- mods_imp_dats[mod_type == "mods_fg" & method == "SMC-FCS Fine-Gray"]$mod[[1]] |>
  cox.zph(terms = TRUE)
prop
plot(prop, resid = FALSE)


prop <- mods_imp_dats[mod_type == "mods_cs1" & method == "SMC-FCS Fine-Gray"]$mod[[1]] |>
  cox.zph(terms = TRUE)
prop
plot(prop, resid = FALSE)
