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

cairo_pdf("analysis/figures/cuminc_mf.pdf", width = 9)
par(family = "Roboto Condensed")
plot(
  np_curves,
  cause = "stacked",
  col = cols[c(4, 6)],
  ylab = "Stacked cumulative incidence",
  lty = c(1, 2),
  xlab = "Time since alloHCT (months)",
  percent = FALSE,
  atrisk.at = seq(0, 60, by = 10),
  atrisk.col = "black",
  atrisk.title = "At risk:",
  legend = FALSE
)
legend(
  x = 0, y = 0.95,
  legend = c("Relapse", "Non-relapse mortality"),
  lty = c(1, 2),
  col = cols[c(4, 6)],
  bty = 'n',
  lwd = c(3, 3)
)
dev.off()
#plot(np_curves, cause = 1, col = "blue") # note the jump in cont prog pats
#plot(np_curves, cause = 2, add = TRUE, lty = 2)


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

library(gtsummary)
tar_load(applied_dat_raw)

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

# Make baseline patient
times_preds <- seq(0.01, 60, by = 0.1)
newpat <- dat[1, ]
newpat[, (sm_predictors) := lapply(.SD, function(col) {
  if (is.factor(col)) levels(col)[1] else 0
}), .SDcols = sm_predictors]

# Compute baseline cumincs + SEs
basecuminc_imp_dats <- mods_imp_dats[mod_type == "mods_fg", .(
  "base" = list({
    obj <- predictCox(
      object = mod[[1]],
      times = times_preds,
      type = "survival",
      newdata = newpat,
      se = TRUE,
      store.iid = "minimal"
    )
    cbind.data.frame(
      "times" = times_preds,
      "cuminc" = 1 - drop(obj$survival),
      "se" = drop(obj$survival.se)
    )
  })
), by = c("tar_batch", "tar_rep", "method")]

# https://rdrr.io/github/survival-lumc/CauseSpecCovarMI/src/R/illustrative-analysis-helpers.R
# Shouldn't the CI be based on t-distribution??
preds_full <- basecuminc_imp_dats[, unlist(base, recursive = FALSE), by = c(
  "tar_batch", "tar_rep", "method"
)]

preds_full[, ':=' (
  p  = cuminc,
  p_trans = cloglog(cuminc),
  var_p = (se)^2
)]
preds_full[, Ui := var_p / (log(1 - p) * (1 - p))^2]
by_vars <- c("method", "times")

preds_summ <- preds_full[, .(
  m = .N,
  Ubar = mean(Ui),
  B = stats::var(p_trans),
  Qbar = mean(p_trans)
), by = by_vars]

preds_summ[, total_var := Ubar + (1 + m^-1) * B, by = by_vars]

preds_final <- preds_summ[, .(
  p_pooled = inv_cloglog(Qbar),
  CI_low = inv_cloglog(Qbar - stats::qnorm(0.975) * sqrt(total_var)),
  CI_upp = inv_cloglog(Qbar + stats::qnorm(0.975) * sqrt(total_var))
), by = by_vars]
preds_final[p_pooled == 0, c("CI_low", "CI_upp") := 0]
preds_final[, CI_width := CI_upp - CI_low]


p_cuminc <- preds_final |>
  ggplot(aes(times, p_pooled, group = method)) +
  #geom_ribbon(aes(ymin = CI_low, ymax = CI_upp, fill = method), alpha = 0.5) +
  geom_step(aes(col = method, linetype = method), linewidth = 0.75) +
  scale_color_manual(values = cols[c(1, 2, 6, 4, 5)]) +
  coord_cartesian(ylim = c(0, 0.225)) +
  labs(
    col = "Method",
    linetype = "Method",
    x = "Time since alloHCT (months)",
    y = "(Pooled) Baseline cumulative incidence"
  )

p_ci_width <- preds_final |>
  ggplot(aes(times, CI_width, group = method)) +
  geom_step(aes(col = method, linetype = method), linewidth = 0.75) +
  scale_color_manual(values = cols[c(1, 2, 6, 4, 5)]) +
  coord_cartesian(ylim = c(0, 0.225)) +
  labs(
    col = "Method",
    linetype = "Method",
    x = "Time since alloHCT (months)",
    y = "(Pooled) 95% Confidence interval width"
  )


p_comb <- p_cuminc + p_ci_width & xlab(NULL) & theme(legend.position = "top")

# p_comb +
#   plot_layout(
#     axis_titles = "collect",
#     guides = "collect"
#   )

p_final <- wrap_elements(panel = p_comb + plot_layout(guides = "collect")) +
  labs(tag = "Time since alloHCT (months)") +
  theme(
    plot.tag = element_text(size = rel(1)),
    plot.tag.position = "bottom"
  )

ggsave(
  plot = p_final,
  filename = "analysis/figures/applied_base_cuminc.pdf",
  width = 11,
  scale = 1,
  height = 7,
  device = cairo_pdf
)


basecuminc_imp_dats[, unlist(base, recursive = FALSE), by = c("tar_batch", "tar_rep", "method")][, .(
  pred = inv_cloglog(mean(cloglog(cuminc)))
), by = c("method", "time")] |>
  ggplot(aes(time, pred, group = method, col = method)) +
  geom_step(aes(linetype = method), linewidth = 1) +
  scale_color_manual(values = cols[c(1, 2, 6, 4, 5)]) +
  coord_cartesian(ylim = c(0, 0.225)) +
  labs(
    col = "Method",
    linetype = "Method",
    x = "Time since alloHCT (months)",
    y = "(Pooled) Baseline cumulative incidence function"
  )

ggsave(
  filename = "analysis/figures/applied_base_cuminc.pdf",
  width = 10,
  scale = 1,
  height = 7,
  device = cairo_pdf
)

# Exponentiate here if needed
pooled_mods <- mods_imp_dats[, .(
  summ = list(tidy(pool(mod), conf.int = TRUE, exponentiate = TRUE))
), by = c("method", "mod_type")][, unlist(summ, recursive = FALSE), by = c("method", "mod_type")]

df_plot <- pooled_mods[!(
  term %in% c("intdiagallo_decades", "year_allo1_decades", "submps_allo1sMF")
) & mod_type == "mods_fg"]

df_plot[, method := factor(
  method,
  levels = c(
    "SMC-FCS Fine-Gray",
    "SMC-FCS cause-spec",
    "MICE cause-spec",
    "MICE subdist",
    "Compl. cases"
  ),
  labels = c(
    "SMC-FCS\nFine-Gray",
    "SMC-FCS\ncause-spec",
    "MICE\ncause-spec",
    "MICE\nsubdist",
    "Complete\ncases"
  )
)]

df_plot[, term := factor(
  term,
  levels = c(
    "ric_allo1reduced",
    "cmv_matchOther",
    "vchromos_preallo1Abnormal",
    "donrel_binOther",
    "hb_allo1",
    "hctci_riskintermediate risk (1-2)",
    "hctci_riskhigh risk (>= 3)",
    "KARNOFSK_threecat80",
    "KARNOFSK_threecat<80",
    "sweat_allo1Yes",
    "age_allo1_decades",
    "PATSEXMale",
    "pb_allo1",
    "ruxo_preallo1yes",
    "tceldepl_binyes",
    "wbc_allo1",
    "WEIGLOSS_allo1Yes"
  ),
  labels = c(
    "Conditioning (reduced)",
    "CMV match (other)",
    "Cytogenetics (abnormal)",
    "Don. relation (other)",
    "Hemoglobin (per 5 g/dL)",
    "HCT-CI (1-2)",
    "HCT-CI (>= 3)",
    "Karnofsky (=80)",
    "Karnofsky (<=70)",
    "Night sweats",
    "Pat. age (decades)",
    "Pat. sex (male)",
    "PB Blasts (per 5%)",
    "Ruxolitinib given",
    "T-cell depletion",
    "WBC count (log)",
    "Weight loss"
  )
)]

df_plot |>
  ggplot(aes(term, estimate, group = method, col = method, shape = method)) +
  geom_point(position = position_dodge(width = 0.75)) +
  geom_pointrange(
    aes(ymin = conf.low, ymax = conf.high),
    position = position_dodge(width = 0.75)
  ) +
  coord_flip() +
  geom_hline(yintercept = 1, linetype = "dotted") +
  scale_color_manual(values = cols[rev(c(1, 2, 6, 4, 5))]) +
  guides(
    colour = ggplot2::guide_legend("Method:", byrow = TRUE, reverse = TRUE),
    shape = ggplot2::guide_legend("Method:", byrow = TRUE, reverse = TRUE)
  ) +
  theme_minimal(base_size = 16, base_family = "Roboto Condensed") +
  theme(legend.position = "top") +
  scale_x_discrete(limits = rev) +
  scale_y_continuous(
   trans = "log",
   breaks = c(0.5, 1, 1.5, 2, 3)
  ) +
  ggstats::geom_stripped_cols(col = NA) +
  theme(
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(hjust = 0),
    axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0))
  ) +
  labs(
    x = "Coefficient",
    y = "Hazard ratio relapse/progression (95% CI)"
  )

 # scale_color_manual(
  #  values = c("")
  #)

ggsave(
  filename = "analysis/figures/applied_forest.pdf",
  width = 8,
  scale = 1.15,
  height = 10,
  device = cairo_pdf
)

modz <- coxph(
  Surv(time_ci_adm, status_ci_adm == 2) ~ sweat_allo1 * ruxo_preallo1,
  data = applied_dat$dat
)

modz

plot(cox.zph(modz))

hi <- copy(applied_dat_raw)
hi[, sweat := forcats::fct_na_value_to_level(sweat_allo1, level = "(missing)")]
hi$sweat
modo <- prodlim(Hist(time_ci_adm, status_ci_adm) ~ sweat,
                data = hi)

df <- copy(applied_dat$dat)
sm_preds <- applied_dat$sm_predictors

colz <- c(
  "ruxo_preallo1",
  "KARNOFSK_threecat",
  "vchromos_preallo1",
  "sweat_allo1",
  "WEIGLOSS_allo1",
  "hctci_risk"
)

df[, (colz) := lapply(.SD, function(col) {
  forcats::fct_na_value_to_level(col, level = "(Missing)")
}), .SDcols = colz]

modo <- prodlim(Hist(time_ci_adm, status_ci_adm) ~ ruxo_preallo1,
                data = df)
plot(modo, cause = 1, legend.cex = 0.75, atrisk = FALSE)


modo <- prodlim(Hist(time_ci_adm, status_ci_adm) ~ KARNOFSK_threecat,
                data = df)
plot(modo, cause = 1, legend.cex = 0.75, atrisk = FALSE)


modo <- prodlim(Hist(time_ci_adm, status_ci_adm) ~ vchromos_preallo1,
                data = df)
plot(modo, cause = 1, legend.cex = 0.75, atrisk = FALSE)


modo <- prodlim(Hist(time_ci_adm, status_ci_adm) ~ sweat_allo1,
                data = df)
plot(modo, cause = 1, legend.cex = 0.75, atrisk = FALSE)

modo <- prodlim(Hist(time_ci_adm, status_ci_adm) ~ WEIGLOSS_allo1,
                data = df)
plot(modo, cause = 1, legend.cex = 0.75, atrisk = FALSE)

modo <- prodlim(Hist(time_ci_adm, status_ci_adm) ~ hctci_risk,
                data = df)
plot(modo, cause = 1, legend.cex = 0.75, atrisk = FALSE)

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

library(kableExtra)
pooled_mods <- mods_imp_dats[, .(
  summ = list(tidy(pool(mod), conf.int = TRUE))
), by = c("method", "mod_type")][, unlist(summ, recursive = FALSE), by = c("method", "mod_type")]

# See table 4 bartlett, make a table for FG1, CSH1 and CSH2 (log HRs)
pooled_mods[, summ := paste0(
  round(estimate, 2), " [",
  round(conf.low, 2), ", ",
  round(conf.high, 2), "]"
)]

pooled_mods[mod_type == "mods_cs1"] |>
  dcast(term ~ method, value.var = "summ") |>
  kbl() |>
  kable_styling("striped")

df_tbl <- dcast(pooled_mods, term + method ~ mod_type, value.var = "summ")
df_tbl[, term := factor(
  term,
  levels = c(
    "ric_allo1reduced",
    "cmv_matchOther",
    "vchromos_preallo1Abnormal",
    "donrel_binOther",
    "hb_allo1",
    "hctci_riskintermediate risk (1-2)",
    "hctci_riskhigh risk (>= 3)",
    "intdiagallo_decades",
    "KARNOFSK_threecat80",
    "KARNOFSK_threecat<80",
    "submps_allo1sMF",
    "sweat_allo1Yes",
    "age_allo1_decades",
    "PATSEXMale",
    "pb_allo1",
    "ruxo_preallo1yes",
    "tceldepl_binyes",
    "wbc_allo1",
    "WEIGLOSS_allo1Yes",
    "year_allo1_decades"
  ),
  labels = c(
    "Conditioning: reduced",
    "CMV match: other",
    "Cytogenetics: abnormal",
    "Don. relation: other",
    "Hemoglobin (per $5$ g/dL)",
    "HCT-CI ($1-2$)",
    "HCT-CI ($\\geq 3$)",
    "Interval diagnosis to alloHCT (decades)",
    "Karnofsky ($=80$)",
    "Karnofsky ($\\leq 70$)",
    "MF subtype: sMF",
    "Night sweats: yes",
    "Pat. age (decades)",
    "Pat. sex: male",
    "PB Blasts (per $5$\\%)",
    "Ruxolitinib given: yes",
    "T-cell depletion: yes",
    "WBC count (log)",
    "Weight loss: yes",
    "Year of alloHCT (decades)"
  )
)]


tb_breaks <- rep(5, length(levels(df_tbl$term)))
names(tb_breaks) <- levels(df_tbl$term)

df_tbl[order(term)][, !c("term")] |>
  kbl(
    format = "latex",
    longtable = TRUE,
    booktabs = TRUE,
    caption = "Longtable",
    align = c("l", "r", "r", "r"),
    col.names = c(
      "Term + method",
      "Relapse subdist.~log HR",
      "Relapse cause-spec.~log HR",
      "NRM cause-spec.~log HR"
    ),
    escape = FALSE
  ) |>
  kable_styling(
    latex_options = c("repeat_header"),
    repeat_header_continued = TRUE,
    repeat_header_method = "replace",
    font_size = 9#,
    #latex_options = c("repeat_header")
  ) |>
  pack_rows(
    index = tb_breaks,
    escape = FALSE
    #extra_latex_after = "\\rowcolor{gray!6}"
    #index =c(" " = 0, "Group 1" = 5, "Group 2" = 0)
  )
  #collapse_rows(columns =1:2, valign = "top")


# Tbl summary
df <- applied_dat$dat
df_predz <- copy(df[, ..sm_predictors]) #|>
# Re-transform for descriptives
df_predz[, ':=' (
  age_allo1_decades = (age_allo1_decades + 6) * 10,
  year_allo1_decades = year_allo1_decades * 10 + 10 + 2009,
  intdiagallo_decades = intdiagallo_decades * 10,
  pb_allo1 = pb_allo1 * 5,
  hb_allo1 = hb_allo1 * 5 + 10,
  wbc_allo1 = exp(wbc_allo1 + log(15.1))
)]

attr(df_predz$hctci_risk, "label") <- "HCT-CI risk category"
attr(df_predz$KARNOFSK_threecat, "label") <- "Karnosfky performance score"
attr(df_predz$tceldepl_bin, "label") <- "T-cell depletion (in- or ev-vivo)"
attr(df_predz$age_allo1_decades, "label") <- "Patient age (years)"
attr(df_predz$cmv_match, "label") <- "Patient/donor CMV match"
attr(df_predz$intdiagallo_decades, "label") <- "Interval diagnosis-transplantation (years)"
attr(df_predz$year_allo1_decades, "label") <- "Year of transplantation"

library(kableExtra)
tbl_summary(
  df_predz,
  missing_text = "(Missing)",
  type = all_dichotomous() ~ "categorical"
) |>
  as_kable_extra(
    format = "latex",
    booktabs = TRUE,
    longtable = TRUE
  )


df <- applied_dat$dat
varz <- c("time_ci_adm", "status_ci_adm", sm_predictors)
df_actual <- df[, ..varz]
df_actual[, CCA_ind := as.numeric(complete.cases(.SD))]
df_actual[, status_ci_adm := factor(status_ci_adm, levels = c(2, 1, 0))]

library(rms)
dd <- datadist(df_actual)
dd$limits$status_ci_adm[2] <- 2
dd$limits$year_allo1_decades[2] <- 0
options(datadist = "dd")

modm <- lrm(#CCA_ind ~ intdiagallo_decades +
  !is.na(KARNOFSK_threecat) ~ intdiagallo_decades +
              rcs(year_allo1_decades, 4) +
              donrel_bin +
              submps_allo1 +
              ric_allo1 +
              cmv_match +
              tceldepl_bin +
              PATSEX +
              age_allo1_decades +
              status_ci_adm *
              rcs(time_ci_adm, 4),
            data = df_actual)

modm
ggplot(Predict(modm, time_ci_adm, status_ci_adm, ref.zero = TRUE))
ggplot(Predict(modm, year_allo1_decades, ref.zero = TRUE))




# Check the CS hazards (make a table) -------------------------------------




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
