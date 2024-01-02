# Possible to-do/ideas:
# - Add proposed smcfcs.finegray to smcfcs() package; add options to kmi settings?
# - Visualise missing mechanism with jitter plots! (+ visualise the settings, together with base cumincs)
# - See {mets} for cifreg() and doubleFGR()


# Manu::get_pal("Hoiho")
cols <- c("#CABEE9", "#7C7189", "#FAE093", "#D04E59", "#BC8E7D", "#2F3D70")


theme_set(
  theme_light(base_size = 16) + #, base_family = "Roboto Condensed") +
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

# Edit these times!
times <- c(seq(0.001, 2.5, length.out = 40),seq(2.6, 10, length.out = 10))
scenarios <- expand.grid(
  mechanism = c("correct_FG", "misspec_FG"),
  p = c(0.15, 0.65),
  stringsAsFactors = FALSE
)
#newdat <- data.frame(X = c(0, 1), Z = c(0, 0))



# Intermediate tests ------------------------------------------------------



compute_true(
  t = times,
  newdat = data.frame(X = 0, Z = 0),
  #params = tar_read(params_weibull_lfps_0.65),
  #model_type = "misspec_FG"#,
  params = tar_read(true_params_correct_FG_0.65),
  model_type = "correct_FG"
) |>
  melt(
    variable.name = "hazard_type",
    value.name = "hazard",
    measure.vars = c("subdist_haz", "cs_haz")
  ) |>
  _[!(hazard_type == "subdist_haz" & cause == 2)] |>
  ggplot(aes(time, hazard)) +
  geom_line(
    aes(
      group = interaction(hazard_type, cause),
      col = interaction(hazard_type, cause),
      linetype = interaction(hazard_type, cause)
    ),
    linewidth = 1.5, alpha = 0.8
  ) +
  scale_color_manual(
    breaks = c("cs_haz.1", "cs_haz.2", "subdist_haz.1"),
    values = cols[1:3],
    labels = c("Cause-spec cause 1", "Cause-spec cause 2", "Subdist. cause 1")
  ) +
  scale_linetype_manual(
    breaks = c("cs_haz.1", "cs_haz.2", "subdist_haz.1"),
    values = 1:3,
    labels = c("Cause-spec cause 1", "Cause-spec cause 2", "Subdist. cause 1")
  ) +
  coord_cartesian(ylim = c(0, 2))







true_vals_ls <- mapply(
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
)

dat_true <- rbindlist(true_vals_ls)
dat_HRs <- dat_true[, .(
  HR_cs = cs_haz[X == 1] / cs_haz[X == 0],
  HR_subdist = subdist_haz[X == 1] / subdist_haz[X == 0]
), by = c("time", "cause", "p", "mech")]

dat_HRs |>
  melt(
    variable.name = "hazard_type",
    value.name = "HR",
    measure.vars = c("HR_cs", "HR_subdist")
  ) |>
  ggplot(aes(time, HR)) +
  geom_line(
    aes(
      group = interaction(hazard_type, cause),
      col = interaction(hazard_type, cause),
      linetype = interaction(hazard_type, cause)
    ),
    linewidth = 1.5, alpha = 0.8
  ) +
  facet_wrap(p * mech ~ ., ncol = 4, labeller = all_labels)


# Start plots
dat_p1 <- dat_true[X == 0 & Z == 0] |>
  melt(
    variable.name = "hazard_type",
    value.name = "hazard",
    measure.vars = c("subdist_haz", "cs_haz")
  )

dat_p1[!(hazard_type == "subdist_haz" & cause == 2)] |>
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
    values = cols[1:3],
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
