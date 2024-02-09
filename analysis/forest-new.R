impdats <- data.table(tar_read(applied_impdats))
tar_load(applied_dat)
options(contrasts = rep("contr.treatment", 2))
impdats[, c("donsex", "patsex") := tstrsplit(sexmatch, " into ")]


impdats[, .(.N), by = c("tar_batch", "tar_rep", "method")]

options(contrasts = rep("contr.treatment", 2))


# Prepare the sms
sm_form_fg <- reformulate(termlabels = c(
  applied_dat$sm_predictors
), response = "Surv(newtimes, newevent)")
sm_form_fg <- update(sm_form_fg, . ~ . - hb_allo1 + I(hb_allo1 / 5) -
                       pb_allo1 + I(10 * pb_allo1 / 5) - log_wbc_allo1 +
                       I(exp(log_wbc_allo1) / 5) - agedonor_allo1_decades
                       )
sm_form_cs1 <- update(sm_form_fg, Surv(time_ci_adm, status_ci_adm == 1) ~ .)
sm_form_cs2 <- update(sm_form_fg, Surv(time_ci_adm, status_ci_adm == 2) ~ .)

applied_dat_raw[!is.na(donrel_bin)] |>
  ggplot(aes(agedonor_allo1_decades, col = donrel)) +
  geom_density()

tar_load(applied_dat_raw)

applied_dat_raw$donrel |> table()

applied_dat_raw[, "donrel_cat" := fcase(
  donrel == "Identical sibling", "Identical sibling",
  donrel == "MMRD", "MMRD",
  default = "Other unrelated"
)]

applied_dat_raw |>
  ggplot(aes(agedonor_allo1_decades, fill = donrel_bin)) +
  geom_density(alpha = 0.5)

coxph(
  Surv(time_ci_adm, status_ci_adm == 1) ~
    donrel_cat +
    agedonor_allo1_decades,
  applied_dat_raw
)

modo <- prodlim(Hist(time_ci_adm, status_ci_adm) ~ donrel_bin, #donrel_cat,
                data = applied_dat_raw[agedonor_allo1_decades > 1.5])

#hist(applied_dat_raw$agedonor_allo1_decades)
plot(modo, cause = 1,
     legend.cex = 0.75,atrisk = FALSE)

coxph(
  Surv(time_ci_adm, status_ci_adm == 1) ~
    age_allo1_decades +
    donrel_cat,
  applied_dat_raw
)

coxph(
  Surv(time_ci_adm, status_ci_adm == 1) ~
    age_allo1_decades +
    donrel_bin,
  applied_dat$dat[agedonor_allo1_decades < 1]
)

coxph(
  Surv(time_ci_adm, status_ci_adm == 1) ~
    age_allo1_decades +
    agedonor_allo1_decades,
  applied_dat$dat[donrel_bin == "Other"]
)

library(gtsummary)
predos <- applied_dat$sm_predictors
applied_dat$dat[, ..predos] |>
  tbl_summary(by = "donrel_bin") |>
  add_p()

test <- impdats[, .(
  mods_fg = list(coxph(sm_form_fg, data = .SD, x = TRUE))#,
  #mods_cs1 = list(coxph(sm_form_cs1, data = .SD, x = TRUE)),
  #mods_cs2 = list(coxph(sm_form_cs2, data = .SD, x = TRUE))
), by = c("tar_batch", "tar_rep", "method")]

test[method == "MICE subdist"][["mods_fg"]] |>
  pool() |>
  howManyImputations::how_many_imputations()

# Fit the models
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

glm(
  is.na(pb_allo1) ~
    age_allo1_decades +
    cmv_match +
    donrel_bin +
    intdiagtr_allo1 +
    year_allo1 +
    time_ci_adm +
    factor(status_ci_adm) +
    ,
  data = applied_dat_raw
) |> summary()

mods_imp_dats[method == "Compl. cases"][["mod"]]

mods_imp_dats

mods_imp_dats[, .(
  summ = list(tidy(pool(mod), conf.int = TRUE))
), by = c("method", "mod_type")][, unlist(
  summ, recursive = FALSE
), by = c("method", "mod_type")] |>
  ggplot(aes(term, estimate, group = method, col = method, shape = method)) +
  geom_point(position = position_dodge(width = 0.75)) +
  geom_pointrange(
    aes(ymin = conf.low, ymax = conf.high),
    position = position_dodge(width = 0.75)
  ) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  facet_wrap(~ mod_type) +
  #facet_grid(mod_type ~ .) +
  scale_color_manual(values = cols[c(1, 2, 6, 4, 5)]) +
  theme(legend.position = "top") +
  scale_x_discrete(limits = rev)# +
  #scale_y_continuous(
  #  trans = "exp",
  #  breaks = log(c(0.5, 0.75, 1, 1.5, 2, 3)),
  #  labels = c(0.5, 0.75, 1, 1.5, 2, 3)
  #)

# Plot subset of variables,
#.. with legend saying model additionally adjust for year of allo,
# intdiagtr, etc.

