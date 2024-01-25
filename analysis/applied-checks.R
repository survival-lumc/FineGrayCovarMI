impdats <- data.table(tar_read(applied_impdats))
tar_load(applied_dat)
options(contrasts = rep("contr.treatment", 2))

#### MAKE THE CS MODELS + FG MODEL AND SAVE IN NESTED DTS??

# For Roboto font and Manu "Hoiho" palette
library(extrafont) # Add to packages file?
cols <- c("#CABEE9", "#7C7189", "#FAE093", "#D04E59", "#BC8E7D", "#2F3D70")
theme_set(
  theme_light(base_size = 16, base_family = "Roboto Condensed") +
    theme(
      strip.background = element_rect(fill = cols[2], colour = "white"),
      strip.text = element_text(colour = 'white')
    )
)

impdats[, .(.N), by = c("id", "method")]
impdats[, .(.N), by = c("tar_batch", "tar_rep", "method")]
# tarbatch * tarrep = imp

# Check cuminc orig data
np_curves <- prodlim(Hist(time_ci_adm, status_ci_adm) ~ 1, data = applied_dat$dat)
plot(np_curves, cause = 1, col = ) # note the jump in cont prog pats
plot(np_curves, cause = 2, add = TRUE, lty = 2)

# Make a stacked plot instead?
library(mstate)
tmat <- trans.comprisk(K = 2)
dat_long <- msprep(
  time = c(NA, "time_ci_adm", "time_ci_adm"),
  status = cbind(
    NA,
    as.numeric(applied_dat$dat$status_ci_adm == 1),
    as.numeric(applied_dat$dat$status_ci_adm == 2)
  ),
  id = "id",
  keep = applied_dat$sm_predictors,
  trans = tmat,
  data = data.frame(applied_dat$dat)
)

dat_long_exp <- expand.covs(
  dat_long,
  covs = applied_dat$sm_predictors,
  longnames = F
)

marg <- coxph(Surv(time, status) ~ strata(trans), data = dat_long_exp)
msf <- msfit(marg, trans = tmat)
plot(msf, use.ggplot = TRUE)
ptrans <- probtrans(msf, predt = 0)
plot(ptrans, use.ggplot = TRUE, ord = c(3, 2, 1), type = "single")

form_mstate <- reformulate(
  response = "Surv(time, status)",
  termlabels = c(
    colnames(dat_long_exp)[20:45],
    colnames(dat_long_exp)[20:45],
    "strata(trans)"
  )
)

coxph(form_mstate, data = dat_long_exp)


sm_form <- reformulate(
  termlabels = applied_dat$sm_predictors,
  response = "Surv(newtimes, newevent)"
)

mods_imp_dats <- impdats[, .(
  mods = list(coxph(sm_form, data = .SD, x = TRUE))
), by = c("tar_batch", "tar_rep", "method")]

times_preds <- seq(0.01, 60, by = 0.1)

hi <- mods_imp_dats[, .(
  "base" = list(
    cbind.data.frame(
      time = times_preds,
      cuminc = 1 - predictCox(object = mods[[1]], times = times_preds, centered = FALSE)$survival
    )
  )
), by = c("tar_batch", "tar_rep", "method")]

hi[, unlist(base, recursive = FALSE), by = c("tar_batch", "tar_rep", "method")][, .(
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


# PAIRWISE PREDS AT 36 months??
horizon <- 36

preds_imp_dats <- impdats[, .(
  id = id,
  preds = {
    mod <- coxph(sm_form, data = .SD, x = TRUE)
    drop(predictRisk(mod, newdata = .SD, times = horizon))
  }
), by = c("tar_batch", "tar_rep", "method")]

df <- applied_dat$dat
id_cc <- df[complete.cases(df), ][["id"]]
preds_cc <- preds_imp_dats[, .(
  preds_pooled = inv_cloglog(mean(cloglog(preds)))
), by = c("method", "id")]#[!(id %in% id_cc)]
preds_cc[id %in% id_cc][method %in% c("Compl. cases", "SMC-FCS Fine-Gray")]
wide <- dcast(
  preds_cc[id %in% id_cc],
  id ~ method, value.var = "preds_pooled"
)

wide |>
  ggplot(aes(`Compl. cases`, `SMC-FCS Fine-Gray`)) +
  geom_point(alpha = 0.5, size = 3) +
  scale_color_viridis_b(direction = -1) +
  geom_abline(a = 0, b = 1, col = "lightblue", linewidth = 1) +
  geom_smooth()

preds_cc[id %in% id_cc] |>
  ggplot(aes(preds_pooled)) +
  geom_density(aes(col = method), alpha = 0.5, fill = NA)

GGally::ggpairs(wide, columns = 2:6)


# COEFS:
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


# Baseline hazards --------------------------------------------------------




# Check imputed cens times ------------------------------------------------



tar_load(applied_dat)
dat_sub <- applied_dat$dat
dat_sub[, year_allo1 := factor(year_allo1)]

cens_imps <- kmi::kmi(
  formula = Surv(time_ci_adm, status_ci_adm != 0) ~ year_allo1 +
    intdiagtr_allo1 +
    age_allo1_decades + PATSEX + ric_allo1 + tbi_allo1 + cmv_match,
  data = data.frame(dat_sub),
  etype = status_ci_adm,
  failcode = 2, # non-relapse mortality is the outcome
  nimp = 50, # make bigger later, set globally
  bootstrap = TRUE,
  nboot = 50
)

imp_dats_cens <- lapply(cens_imps$imputed.data, function(imp_dat) {
  cbind.data.frame(cens_imps$original.data, imp_dat)
}) |> rbindlist(idcol = ".imp")

imp_dats_cens[, H1_subdist := compute_marginal_cumhaz(
  timevar = newtimes,
  statusvar = newevent,
  cause = 2,
  type = "cause_spec"
), by = .imp]
imp_dats_cens[, .imp := factor(.imp)]

imp_dats_cens |>
  ggplot(aes(newtimes, H1_subdist)) +
  geom_step(aes(group = .imp, col = .imp)) +
  theme(legend.position = "none")



# Try experiment with simulated data --------------------------------------


dat <- generate_dataset(
  n = 1e3,
  list(
    mechanism = "correct_FG",
    params = tar_read(true_params_correct_FG_0.65),
    censoring_type = "exponential",
    censoring_params = list("exponential" = "0.49")
    #censoring_params = list("exponential" = "0.49 * exp(Z)")
  ),
  list(mech_params = list("prob_missing" = 0, "mechanism_expr" = "1"))
)

table(dat$D)
#coxph(Surv(time, D == 0) ~ Z, data = dat)

cens_imps <- kmi::kmi(
  formula = Surv(time, D != 0) ~ 1, # ~ Z
  data = data.frame(dat),
  etype = D,
  failcode = 1, # non-relapse mortality is the outcome
  nimp = 50, # make bigger later, set globally
  bootstrap = F,
  nboot = 100
)

imp_dats_cens <- lapply(cens_imps$imputed.data, function(imp_dat) {
  cbind.data.frame(cens_imps$original.data, imp_dat)
}) |> rbindlist(idcol = ".imp")

imp_dats_cens[, H1_subdist := compute_marginal_cumhaz(
  timevar = newtimes,
  statusvar = newevent,
  cause = 1,
  type = "cause_spec"
), by = .imp]
imp_dats_cens[, .imp := factor(.imp)]

dat[, H1_subdist := compute_marginal_cumhaz(
  timevar = time,
  statusvar = D,
  cause = 1,
  type = "subdist"
)]

dat[, "V_true" := (D == 1) * time + (D != 1) * cens_time]
dat

dat[, H1_subdist_true := compute_marginal_cumhaz(
  timevar = V_true,
  statusvar = D,
  cause = 1,
  type = "cause_spec"
)]

p_subdens <- dat[D == 1] |>
  ggplot(aes(time)) +
  geom_histogram(fill = "lightblue", col = "black", bins = 50) +
  coord_cartesian(xlim = c(0, 10)) +
  theme_minimal()

p_subdist <- imp_dats_cens |>
  ggplot(aes(newtimes, H1_subdist)) +
  geom_step(aes(group = .imp), alpha = 0.75, col = "gray") +
  theme(legend.position = "none") +
  #coord_cartesian(ylim = c(0, 0.5), xlim = c(0, 10)) +
  # Add the weighted one
  geom_step(data = dat, aes(x = time), linewidth = 1.5, linetype = "dashed") +
  geom_step(
    data = dat,
    aes(x = V_true, y = H1_subdist_true),
    linewidth = 1.5,
    col = "red",
    linetype = "dotted"
  ) +
  theme_light()

p_subdist / p_subdens + plot_layout(heights = c(3, 1))
