# Possible to-do/ideas:
# - Add proposed smcfcs.finegray to smcfcs() package; add options to kmi settings?
# - Visualise missing mechanism with jitter plots! (+ visualise the settings, together with base cumincs)
# - See {mets} for cifreg() and doubleFGR()


# Manu::get_pal("Hoiho")
# need to load extrafont if we want Roboto Condensed
library(extrafont)
cols <- c("#CABEE9", "#7C7189", "#FAE093", "#D04E59", "#BC8E7D", "#2F3D70")

theme_set(
  theme_light(base_size = 16, base_family = "Roboto Condensed") +
    theme(
      strip.background = element_rect(fill = Manu::get_pal("Hoiho")[[2]], colour = "white"),
      strip.text = element_text(colour = 'white')
    )
)

# Global labels
pspace_labels <- c("0.15" = "p = 0.15", "0.65" = "p = 0.65")
failure_model_labels <- c("correct_FG" = "Well specified FG", "misspec_FG" = "Misspecified FG")
all_labels <- labeller(p = pspace_labels, mech = failure_model_labels)


# Figure 2 ----------------------------------------------------------------


times <- c(
  seq(0.001, 2.5, length.out = 400),
  seq(2.6, 10, length.out = 100)
)
scenarios <- expand.grid(
  mechanism = c("correct_FG", "misspec_FG"),
  p = c(0.15, 0.65),
  stringsAsFactors = FALSE
)
newdat <- data.frame(X = c(0, 1), Z = c(0, 0))

dat_true <- mapply(
  function(mechanism, p) {
    params_slug <- switch(
      mechanism,
      "misspec_FG" = paste0("params_weibull_lfps_", p),
      "correct_FG" = paste0("true_params_correct_FG_", p)
    )
    params <- do.call(tar_read, list(as.symbol(params_slug)))
    compute_true(
      t = times,
      newdat = newdat,
      params = params,
      model_type = mechanism
    ) |> cbind(p = p, mech = mechanism)
  },
  mechanism = scenarios$mechanism,
  p = scenarios$p,
  SIMPLIFY = FALSE,
  USE.NAMES = FALSE
) |> rbindlist()



# Start plots
dat_p1 <- dat_true[X == 0 & Z == 0] |>
  melt(
    variable.name = "hazard_type",
    value.name = "hazard",
    measure.vars = c("subdist_haz", "cs_haz")
  )

p1 <- dat_p1[!(hazard_type == "subdist_haz" & cause == 2)] |>
  ggplot(aes(time, hazard)) +
  geom_line(
    aes(
      group = interaction(hazard_type, cause),
      col = interaction(hazard_type, cause),
      linetype = interaction(hazard_type, cause)
    ),
    linewidth = 1.5, alpha = 0.8
  ) +
  facet_wrap(p * mech ~ ., ncol = 4, labeller = all_labels) +
  scale_color_manual(
    breaks = c("cs_haz.1", "cs_haz.2", "subdist_haz.1"),
    values = cols[c(1, 2, 5)],
    labels = c("Cause-spec cause 1", "Cause-spec cause 2", "Subdist. cause 1")
  ) +
  scale_linetype_manual(
    breaks = c("cs_haz.1", "cs_haz.2", "subdist_haz.1"),
    values = 1:3,
    labels = c("Cause-spec cause 1", "Cause-spec cause 2", "Subdist. cause 1")
  ) +
  coord_cartesian(ylim = c(0, 2)) +
  labs(
    y = "Baseline hazards",
    x = "Time",
    linetype = NULL,
    col = NULL
  )

p2 <- dat_true[X == 0 & Z == 0] |>
  ggplot(aes(time, cuminc)) +
  geom_line(
    aes(group = cause, col = cause, linetype = cause),
    linewidth = 1.5, alpha = 0.8
  ) +
  facet_wrap(p * mech ~ ., ncol = 4) +
  scale_color_manual(
    values = cols[c(1, 2)],
    labels = c("Cause 1", "Cause 2")
  ) +
  scale_linetype_manual(
    values = 1:2,
    labels = c("Cause 1", "Cause 2")
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) +
  labs(
    y = "Baseline cumulative incidence",
    x = "Time",
    linetype = NULL,
    col = NULL
  )

# Lower panel
dat_p3 <- dat_true[, .(
  HR_cs = cs_haz[X == 1] / cs_haz[X == 0],
  HR_subdist = subdist_haz[X == 1] / subdist_haz[X == 0]
), by = c("time", "cause", "p", "mech")] |>
  melt(
    variable.name = "HR_type",
    value.name = "HR",
    measure.vars = c("HR_cs", "HR_subdist")
  )

p3 <- dat_p3[!(HR_type == "HR_subdist" & cause == 2)]|>
  ggplot(aes(time, log(HR))) +
  geom_line(
    aes(
      group = interaction(HR_type, cause),
      col = interaction(HR_type, cause),
      linetype = interaction(HR_type, cause)
    ),
    linewidth = 1.5,
    alpha = 0.8
  ) +
  facet_wrap(p * mech ~ ., ncol = 4) +
  scale_color_manual(
    breaks = c("HR_cs.1", "HR_cs.2", "HR_subdist.1"),
    values = cols[c(1, 2, 5)],
    labels = c("Cause-spec cause 1", "Cause-spec cause 2", "Subdist. cause 1")
  ) +
  scale_linetype_manual(
    breaks = c("HR_cs.1", "HR_cs.2", "HR_subdist.1"),
    values = 1:3,
    labels = c("Cause-spec cause 1", "Cause-spec cause 2", "Subdist. cause 1")
  ) +
  geom_hline(yintercept = 0, col = "black", linetype = "dotted") +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) +
  labs(
    y = "Log hazard ratio X = 1 vs X = 0",
    x = "Time",
    linetype = NULL,
    col = NULL
  )

p1 / p2 / p3

ggsave(
  plot = p1 / p2 / p3,
  "analysis/figures/scenarios_vis.pdf", # pdf later
  dpi = 300,
  width = 13,
  height = 9,
  device = cairo_pdf
)




# Other tests -------------------------------------------------------------


p_var <- seq(0.05, 0.95, by = 0.1)

df_ls <- lapply(p_var, function(p) {
  test <- compute_true(
    t = times,
    newdat = newdat,
    params = list(
      "cause1" = list(
        "formula" = ~ X + Z,
        "betas" = c(0.75, 0.5),
        "p" = p,
        "base_rate" = 1,
        "base_shape" = 0.75
      ),
      "cause2" = list(
        "formula" = ~ X + Z,
        "betas" = c(0.75, 0.5),
        "base_rate" = 1,
        "base_shape" = 0.75
      )
    ),
    model_type = "correct_FG"
  )

  test[, .(
    HR_cs = cs_haz[X == 1] / cs_haz[X == 0],
    HR_subdist = subdist_haz[X == 1] / subdist_haz[X == 0]
  ), by = c("time", "cause")] |>
    melt(
      variable.name = "HR_type",
      value.name = "HR",
      measure.vars = c("HR_cs", "HR_subdist")
    ) |>
    cbind("p" = p)
})

rbindlist(df_ls) |>
  ggplot(aes(time, HR)) +
  geom_line(
    aes(
      group = interaction(HR_type, cause),
      col = interaction(HR_type, cause),
      linetype = interaction(HR_type, cause)
    ),
    linewidth = 1.5,
    alpha = 0.8
  ) +
  facet_grid(cause ~ p, scales = "free")

# Second attempt


df_ls_misspec <- lapply(p_var, function(p) {
  pars <- list(
    "cause1" = list(
      "formula" = ~ X + Z,
      "betas" = c(0.75, 0.5),
      "p" = p,
      "base_rate" = 1,
      "base_shape" = 0.75
    ),
    "cause2" = list(
      "formula" = ~ X + Z,
      "betas" = c(0.75, 0.5),
      "base_rate" = 1,
      "base_shape" = 0.75
    )
  )

  pars_misspec <- recover_weibull_lfps(
    large_dat = generate_dataset(
      n = 5000,
      args_event_times = list(
        mechanism = "correct_FG",
        params = pars,
        censoring_type = "none"
      ),
      args_missingness = list(mech_params = list("prob_missing" = 0))
    ),
    params_correct_FG = pars
  )

  test <- compute_true(
    t = times,
    newdat = newdat,
    params = pars_misspec,
    model_type = "misspec_FG"
  )

  test[, .(
    HR_cs = cs_haz[X == 1] / cs_haz[X == 0],
    HR_subdist = subdist_haz[X == 1] / subdist_haz[X == 0]
  ), by = c("time", "cause")] |>
    melt(
      variable.name = "HR_type",
      value.name = "HR",
      measure.vars = c("HR_cs", "HR_subdist")
    ) |>
    cbind("p" = p)
})

rbindlist(df_ls_misspec) |>
  ggplot(aes(time, HR)) +
  geom_line(
    aes(
      group = interaction(HR_type, cause),
      col = interaction(HR_type, cause),
      linetype = interaction(HR_type, cause)
    ),
    linewidth = 1.5,
    alpha = 0.8
  ) +
  facet_grid(cause ~ p, scales = "free")



# Smooth schoenfeld -------------------------------------------------------

pars <- list(
  "cause1" = list(
    "formula" = ~ X + Z,
    "betas" = c(0.75, 0.5),
    "p" = 0.15,
    "base_rate" = 1,
    "base_shape" = 0.75
  ),
  "cause2" = list(
    "formula" = ~ X + Z,
    "betas" = c(0.75, 0.5),
    "base_rate" = 1,
    "base_shape" = 0.75
  )
)

df <- generate_dataset(
  n = 50000,
  args_event_times = list(
    mechanism = "correct_FG",
    params = pars,
    censoring_type = "none"
  ),
  args_missingness = list(mech_params = list("prob_missing" = 0))
)

test <- compute_true(
  t = times,
  newdat = newdat,
  params = pars,
  model_type = "correct_FG"
)

test_long <- test[, .(
  HR_cs = cs_haz[X == 1] / cs_haz[X == 0],
  HR_subdist = subdist_haz[X == 1] / subdist_haz[X == 0]
), by = c("time", "cause")] |>
  melt(
    variable.name = "HR_type",
    value.name = "HR",
    measure.vars = c("HR_cs", "HR_subdist")
  )


with(
  test_long[cause == 2 & HR_type == "HR_cs"],
  {
    plot(time, log(HR), type = "l", lwd = 2)#, ylim = c(0.5, 3))
    abline(a = 1, b = 0, lty = 2)
  }
)

mod <- coxph(Surv(time, D == 2) ~ X + Z, data = df)
zph <- cox.zph(mod, terms = FALSE, transform = "identity")
cox.zph
plot(zph, var = "X1", ylim = c(-6, 6), resid = T, df = 3)
with(
  test_long[cause == 2 & HR_type == "HR_cs"],
  {
    lines(time, log(HR), type = "l", lwd = 2, ylim = c(0.5, 2))
    abline(a = 1, b = 0, lty = 3)
  }
)

data.frame("time" = zph$time, "resid" = zph$y[, 1]) |>
  data.table() |>
  _[time < 10] |>
  ggplot(aes(time, resid)) +
  geom_point(alpha = 0.1) +
  coord_cartesian(ylim = c(-6, 6)) +
  geom_smooth(
    #formula = y ~ splines::ns(x, df = 3),
    #method = "loess",
    #span = 0.9
  ) +
  geom_line(
    data = test_long[cause == 2 & HR_type == "HR_cs"],
    aes(y = log(HR)),
    linewidth = 1
  )



data.frame("time" = zph$time, "resid" = zph$y[, 1]) |>
  ggplot(aes(time, resid)) +
  #geom_point() +
  geom_smooth() +
  geom_line(
    data = test_long[cause == 1 & HR_type == "HR_cs"],
    aes(y = log(HR))
  ) +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 6))

plot(zph$time, zph$y[, 1])
plot(zph, hr = FALSE, var = "X1", resid = FALSE)#, xlim = c(0, 10))
