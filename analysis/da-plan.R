source(here::here("packages.R"))
invisible(lapply(list.files(here::here("R"), full.names = TRUE), source))
tar_load(applied_dat_raw)

hist(qlogis((applied_dat_raw$pb_allo1 + 0.1) / 100))


# Overall missing vals
naniar::gg_miss_var(applied_dat_raw, show_pct = TRUE)
print(naniar::miss_var_summary(applied_dat_raw), n = 60)

#
predictors_def <- c(
  "intdiagtr_allo1", # divide by 10 (so in decades)
  "year_allo1", # divide by 10 (in decades)
  "tbi_allo1",
  "ric_allo1",
  "PATSEX", # maybs replace by sexmatch
  "age_allo1_decades",
  "submps_allo1", # good auxiliary
  "WEIGLOSS_allo1",
  "donrel_bin",
  "cmv_match",
  "KARNOFSK_threecat",
  "wbc_allo1", # log(wbc + 0.1)
  "hb_allo1",
  "vchromos_preallo1",
  "hctci_risk" # since we want variables predictive of NRM too
)

# eventually remove tbi.. since only 7% yes

# But probably not
optional <- c(
  "splenec", # (+ other spleen variables, since too many missings), and not easy to impute..
  # .. because NA for spleen size (should be zero) when splenectomized
  "ruxo_preallo1", # pre-Tx intervention, unsure if only from certain year, this one is possible
  "pb_allo1" # super skewed, could show for misspecification of X | Z?
)

# Optional to (probablt) add:
# - Pat/Don sex match
# - CMV match
# - TBI x RIC combo
# - source!

# Not possible/avoid for now (without asking Linda):
# - Donor age (would need to contact Linda) since we already have continuous (disease) vars to impute


coxph(
  reformulate(
    termlabels = predictors_def,
    response = "Surv(time_ci_adm, status_ci_adm == 1)"
  ),
  data = applied_dat_raw
)


library(gtsummary)
applied_dat_raw[, !c("id_new", "CENTRE_allo1")] |>
  tbl_summary(by = "splenec") |>
  add_p()


tbl_summary(applied_dat_raw, by = "splenec")

cand_preds <- c(
  "submps_allo1",
  "WEIGLOSS_allo1",
  "splenec",
  #"wbc_allo1",
  #"hb_allo1",
  "vchromos_preallo1",
  "ruxo_preallo1",
  "hctci_risk",
  "KARNOFSK_threecat",
  "ric_allo1",
  "tbi_allo1",
  "donrel_bin",
  "age_threecat",
  #"year_allo1",
  #"intdiagtr_allo1",
  "PATSEX",
  "cmv_match"
)

glm(
  is.na(spleen_physicalV5) ~ PATSEX +
    cmv_match +
    year_allo1 +
    intdiagtr_allo1 +
    age_allo1_decades +
    submps_allo1 +
    time_ci_adm +
    factor(status_ci_adm) +
    donrel_bin +
    tbi_allo1 +
    ric_allo1 +
    KARNOFSK_threecat, data = applied_dat_raw
) |>
  vcov() |>
  cov2cor() |>
  round(digits = 3) |>
  corrplot::corrplot()

for (i in seq_along(cand_preds)) {
  form <- reformulate(
    termlabels = cand_preds[i],
    response = "Hist(time_ci_adm, status_ci_adm)"
  )
  cpr <- cmprsk::cuminc(
    applied_dat_raw$time_ci_adm,
    applied_dat_raw$status_ci_adm,
    applied_dat_raw[[cand_preds[i]]]
  )
  mod <- prodlim(form, data = applied_dat_raw)
  plot(
    mod,
    cause = 1,
    atrisk = FALSE,
    lty = c(1, 2),
    legend.cex = 0.9#,
    #atRisk.cex = 0.9
  )
  title(paste0("p = ", round(cpr$Tests[1, "pv"], 4)))
}

sf <- survfit(Surv(time_ci_adm, factor(status_ci_adm)) ~ PATSEX, data = applied_dat_raw)


table(applied_dat_raw$splenec)

prop.table(table(is.na(applied_dat_raw$status_ci_cgvh)))

table(applied_dat_raw$status_ci_agvh2to4)

names(applied_dat_raw)


