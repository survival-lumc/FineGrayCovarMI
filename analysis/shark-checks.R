#tar_meta(fields = c(name, time, bytes, seconds, warnings, error)) |> View()

# (!) Make this into rmd so that all figures are on the github!


# Overview main sims ------------------------------------------------------


# For Roboto font and Manu "Hoiho" palette
library(extrafont) # Add to packages file?
cols <- c("#CABEE9", "#7C7189", "#FAE093", "#D04E59", "#BC8E7D", "#2F3D70")

# Set global ggplot settings
theme_set(
  theme_light(base_size = 16, base_family = "Roboto Condensed") +
    theme(
      strip.background = element_rect(fill = cols[2], colour = "white"),
      strip.text = element_text(colour = 'white')
    )
)

# Global labels for facets
pspace_labels <- c("0.15" = "p = 0.15", "0.65" = "p = 0.65")
failure_model_labels <- c("correct_FG" = "Well specified FG", "misspec_FG" = "Misspecified FG")
censoring_labels <- c(
  "admin" = "Admin. censoring",
  "exponential" = "Random censoring",
  "none" = "No censoring"
)
all_labels <- labeller(
  prob_space = pspace_labels,
  censoring_type = censoring_labels,
  failure_time_model = failure_model_labels
)

# Unnest the whole damn thing
reps_per_scen <- tar_read(simulations_main)[, .(
  coefs = list(rbindlist(coefs_summary, idcol = "rep_id")),
  preds = list(rbindlist(preds_summary, idcol = "rep_id"))
), by = c("prob_space","failure_time_model", "censoring_type", "method")]

# Add proper labels
reps_per_scen[, ':=' (
  censoring_type = factor(censoring_type, levels = c("none", "exponential", "admin")),
  method = factor(
    method,
    levels = c("full", "CCA", "mice_comp", "mice_subdist", "smcfcs_comp", "smcfcs_finegray"),
    labels = c(
      "Full data",
      "Compl. cases",
      "MI cause-spec",
      "MI subdist",
      "SMC-FCS cause-spec",
      "SMC-FCS Fine-Gray"
    )
  )
)]

# Df just for regression coefficients
coefs_main <- rbindlist(
  with(
    reps_per_scen,
    Map(
      f = cbind,
      method = method,
      coefs,
      prob_space = prob_space,
      failure_time_model = failure_time_model,
      censoring_type = censoring_type
    )
  ), fill = TRUE
)

# Some checks
coefs_main[, .(.N), by = c(
  "prob_space",
  "failure_time_model",
  "censoring_type",
  "method",
  "term"
)]

# Do some cleaning
coefs_main[, term := ifelse(grepl(pattern = "^X", term), "X", as.character(term))]


# Performance measures in one go ------------------------------------------


sim_summ <- rsimsum::multisimsum(
  data = coefs_main,
  par = "term",
  estvarname = "estimate",
  se = "std.error",
  true = "true",
  methodvar = "method",
  ref = "Full data",
  ci.limits = c("conf.low", "conf.high"),
  x = TRUE,
  by = c("censoring_type", "failure_time_model", "prob_space")
)


# Jitter plot relative bias -----------------------------------------------


coefs_main[term == "X" & !(method %in% c( "Full data"))] |>
  ggplot(aes(method, 100 * (estimate - true) / true)) +
  geom_jitter(aes(col = method), size = 2.5, width = 0.25, alpha = 0.25, shape = 16) +
  facet_grid(
    failure_time_model * censoring_type ~ prob_space,
    labeller = all_labels
  ) +
  coord_flip() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "none"
  ) +
  scale_y_continuous(minor_breaks = NULL, breaks = c(0, 10, 25, 50, -10, -25, -50)) +
  scale_x_discrete(limits = rev) +
  geom_hline(
    aes(yintercept = 0),
    linetype = "solid",
    linewidth = 1,
    col = cols[3],
    alpha = 1
  ) +
  stat_summary(
    fun = mean,
    geom = "crossbar",
    fatten = 1.5,
    #linewidth = 0.5,
    col = "black",
    alpha = 0.75,
    lineend = "round"
  ) +
  scale_color_manual(values = cols[c(1, 2, 6, 4, 5)]) +
  labs(x = "Method", y = "100 * (Estimate - True) / True (%)")

ggsave(
  filename = "analysis/figures/bias_X.pdf",
  width = 7,
  scale = 1.25,
  height = 10,
  device = cairo_pdf
)


# X variance/coverage summary ---------------------------------------------


crit <- qnorm((1 - 0.95) / 2, lower.tail = FALSE)
df_summ <- data.table(sim_summ$summ)
df_summ_var <- df_summ[stat %in% c("modelse", "empse", "cover")]
df_summ_var[, stat := factor(
  stat, levels = c("empse", "modelse", "cover"),
  labels = c("Emp. SE", "Mod. SE", "Coverage")
)]
df_summ_var[, stat_lab := paste0(
  "bold('", stat, ":')~",
  round(est, 2), "~(", ifelse(mcse < 0.001, "~'<'~0.001", round(mcse, 3)) , ")"
)]



# Bias --------------------------------------------------------------------


# Preb labels
df_lab_bias <- data.table(sim_summ$summ)[stat %in% c("modelse", "empse")]
df_lab_bias[, stat := factor(stat, levels = c("empse", "modelse"), labels = c("Emp. SE", "Mod. SE"))]
df_lab_bias[, stat_lab := paste0(
  "bold('", stat, ":')~",
  round(est, 2), "~(", ifelse(mcse < 0.001, "~'<'~0.001", round(mcse, 3)) , ")"
)]

test <- df_lab_bias[, .(
  lab = paste0("atop(", paste(stat_lab, collapse = ","), ")")
  #lab = paste0("atop(")
), by = c(
  "method",
  "censoring_type",
  "failure_time_model",
  "prob_space",
  "term"
)]




# Master lolly ------------------------------------------------------------

df_lolly <- data.table(sim_summ$summ)[
  term == "X" & !(method %in% c("Full data", "Compl. cases"))
]
dodge_w <- 0.75

left_lab <- df_lolly[stat %in% c("modelse", "empse")][, .(
  left_lab = max(est)
)]#, by = prob_space]

#
df_lab <- df_lolly[stat == "cover"]
df_lab[, stat_lab := paste0(
  "bold('Coverage:')~", round(est, 3), "~(", round(mcse, 4), ")"
)]


df_lolly[stat %in% c("modelse", "empse")] |>
  ggplot(aes(method, est, col = method)) +
  geom_text(
    data = cbind(df_lab, left_lab),#merge(df_lab, left_lab),
    aes(label = stat_lab, y = left_lab + 0.02),
    hjust = 0,
    parse = TRUE,
    family = "Roboto Condensed",
    size = 3
  ) +
  geom_linerange(
    aes(xmin = method, xmax = method, ymin = 0, ymax = est, linetype = stat),
    position = position_dodge(width = dodge_w)
  ) +
  geom_point(
    size = 1.5,
    aes(shape = stat),
    position = position_dodge(width = dodge_w)
  ) +
  geom_point(
    aes(y = est + crit * mcse, group = stat),
    shape = 41,
    size = 1.5,
    position = position_dodge(width = dodge_w)
  ) +
  geom_point(
    aes(y = est - crit * mcse, group = stat),
    shape = 40,
    size = 1.5,
    position = position_dodge(width = dodge_w)
  ) +
  scale_x_discrete(limits = rev) +
  facet_grid(
    failure_time_model * censoring_type ~ prob_space,
    labeller = all_labels
  ) +
  coord_flip(ylim = c(0.05, 0.325)) +
  guides(colour = "none") +
  theme(legend.position = "top") +
  scale_shape_manual(
    "Standard error:",
    values = c(16, 17),
    labels = c("Empirical", "Model-based")
  ) +
  scale_linetype_manual(
    "Standard error:",
    values = c(1:2),
    labels = c("Empirical", "Model-based")
  ) +
  scale_color_manual(values = Manu::get_pal("Hoiho")[c(1, 2, 6, 4)])

ggsave(
  filename = "analysis/figures/perf_X.pdf",
  width = 7,
  scale = 1.25,
  height = 10,
  device = cairo_pdf
)


#lolly_coverage <-
df_lolly[stat == "cover"] |>
  ggplot(aes(method, est, col = method)) +
  geom_linerange(aes(xmin = method, xmax = method, ymin = 0.6, ymax = est)) +
  geom_point(size = 3) +
  geom_point(
    aes(y = est + crit * mcse),
    shape = 41,
    size = 2.5
  ) +
  geom_point(
    aes(y = est - crit * mcse),
    shape = 40,
    size = 2.5
  ) +
  scale_x_discrete(limits = rev) +
  facet_grid(
     failure_time_model * censoring_type * prob_space ~ stat,
    labeller = all_labels,
    scales = "free"
    #scales = "free_y"
  )  +
  theme(
    #axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    panel.grid.major.y = element_blank()
  ) +
  scale_color_manual(values = Manu::get_pal("Hoiho")) +
  geom_hline(yintercept = 0.95, linetype = "dotted") +
  coord_flip(ylim = c(0.65, 1.025)) +
  labs(y = "Coverage (95% CI with MCSE)")


# Coverage ----------------------------------------------------------------


crit <- qnorm((1 - 0.95) / 2, lower.tail = FALSE)

# Coverage
data.table(sim_summ$summ)[stat == "cover" & term == "X"] |>
  ggplot(aes(method, est, col = method)) +
  geom_linerange(aes(xmin = method, xmax = method, ymin = 0.6, ymax = est)) +
  geom_point(size = 3) +
  geom_point(
    aes(y = est + crit * mcse),
    shape = 41,
    size = 2.5
  ) +
  geom_point(
    aes(y = est - crit * mcse),
    shape = 40,
    size = 2.5
  ) +
  scale_x_discrete(limits = rev) +
  facet_grid(
    prob_space ~ failure_time_model * censoring_type,
    labeller = all_labels,
    scales = "free_y"
  )  +
  theme(
    #axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    panel.grid.major.y = element_blank()
  ) +
  scale_color_manual(values = Manu::get_pal("Hoiho")) +
  geom_hline(yintercept = 0.95, linetype = "dotted") +
  coord_flip(ylim = c(0.65, 1.025)) +
  labs(y = "Coverage (95% CI with MCSE)")


# Zippers -----------------------------------------------------------------


coefs_bis <- copy(coefs_main)
coefs_bis[, ':=' (
  estimate = estimate - true,
  conf.low = conf.low - true,
  conf.high = conf.high - true
)]

sim_summ_zip <- simsum(
  data = coefs_bis[term == "X"],
  estvarname = "estimate",
  se = "std.error",
  true = 0,
  methodvar = "method",
  ref = "Full data",
  ci.limits = c("conf.low", "conf.high"),
  x = TRUE,
  by = c("censoring_type", "failure_time_model", "prob_space")
)

autoplot(sim_summ_zip, type = "zip")


# Try the predictions now -------------------------------------------------


preds_main <- rbindlist(
  with(
    reps_per_scen,
    Map(
      f = cbind,
      method = method,
      preds,
      prob_space = prob_space,
      failure_time_model = failure_time_model,
      censoring_type = censoring_type
    )
  ), fill = TRUE
)
preds_main[is.na(imp), imp := 0]
new_pat <- c(1, 1) #X is only 0, 1

preds_main[, "lp" := sum(unlist(coefs) * new_pat), by = c(
  "method",
  "imp",
  "time",
  "rep_id",
  "prob_space",
  "failure_time_model",
  "censoring_type"
)]

preds_main[, pred := 1 - (1 - base_cuminc)^exp(lp)]

pooled_preds <- preds_main[, .(pooled_pred = inv_cloglog(mean(cloglog(pred)))), by = c(
  "method",
  "rep_id",
  "time",
  "prob_space",
  "failure_time_model",
  "censoring_type"
)]
pooled_preds[, ':=' (X = 1, Z = 1)]

true_vals <- tar_read(true_cuminc_all)[cause == 1]
preds_df <- merge(pooled_preds, true_vals)

sim_summ_preds <- rsimsum::simsum(
  data = preds_df,
  estvarname = "pooled_pred",
  se = NULL,
  true = "cuminc",
  methodvar = "method",
  ref = "Full data",
 # x = TRUE,
  by = c("censoring_type", "failure_time_model", "prob_space", "time")
)

#preds_df[is.na(pooled_pred)]
summ_df <- data.table(sim_summ_preds$summ)
summ_df[, time := as.numeric(as.character(time))]
summ_df$stat |> unique()

summ_df[stat == "thetamean" & !(method %in% "Full data")] |>
  ggplot(aes(time, 100 * est)) +
  geom_line(aes(col = method, group = method, linetype = method), size = 2, alpha = 0.7) +
  #geom_line(aes(y = 100 * cuminc), col = "black", size = 1.25) +
  facet_grid(
    prob_space ~ failure_time_model * censoring_type,
    scales = "free",
    labeller = all_labels
  ) +
  #coord_cartesian(ylim = c(-5, 5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = Manu::get_pal("Hoiho")) +
  labs(x = "Time", y = "Bias baseline cumulative incidence (%)")


summ_df[stat == "mse" & !(method %in% "Full data")] |>
  ggplot(aes(time, 100 * sqrt(est))) +
  geom_line(aes(col = method, group = method, linetype = method), size = 2, alpha = 0.7) +
  #geom_line(aes(y = 100 * cuminc), col = "black", size = 1.25) +
  facet_grid(
    prob_space ~ failure_time_model * censoring_type,
    scales = "free",
    labeller = all_labels
  ) +
  #coord_cartesian(ylim = c(-5, 5)) +
  #geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = Manu::get_pal("Hoiho")) +
  labs(x = "Time", y = "RMSE (%)")

summ_df[stat == "empse" & !(method %in% "Full data")] |>
  ggplot(aes(time, est)) +
  geom_line(aes(col = method, group = method, linetype = method), size = 2, alpha = 0.7) +
  #geom_line(aes(y = 100 * cuminc), col = "black", size = 1.25) +
  facet_grid(
    prob_space ~ failure_time_model * censoring_type,
    scales = "free",
    labeller = all_labels
  ) +
  #coord_cartesian(ylim = c(-5, 5)) +
  #geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = Manu::get_pal("Hoiho")) +
  labs(x = "Time", y = "RMSE (%)")



# Try the pairwise tings
preds_df[method %in% c("SMC-FCS cause-spec", "SMC-FCS Fine-Gray")] |>
  dcast(
    formula = rep_id + cuminc + prob_space + failure_time_model + censoring_type ~ method,
    value.var = "pooled_pred"
  ) |>
  ggplot(aes(`SMC-FCS cause-spec`, `SMC-FCS Fine-Gray`)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  facet_grid(
    prob_space ~ failure_time_model * censoring_type,
    scales = "fixed"
  ) +
  geom_smooth()



# What are the true ones? -------------------------------------------------


preds_df[censoring_type != "admin"][, .(
  .N,
  pred = mean(pooled_pred, na.rm = TRUE)
), by = c(
  "time",
  "method",
  "prob_space",
  "failure_time_model",
  "cuminc",
  "censoring_type",
  "X",
  "Z"
)] |>  #[!(method %in% c("full", "mice_comp"))] |>
  ggplot(aes(time, pred)) +
  geom_line(aes(linetype = method, col = method, group = method), linewidth = 1) +
  geom_line(aes(y = cuminc), col = "black", linewidth = 1) +
  facet_grid(
    failure_time_model * censoring_type ~ prob_space,
    labeller = all_labels
  )

#[, .(


preds_df[censoring_type != "admin" & method != "full"][, .(
  .N,
  rmse = sqrt(mean((pooled_pred - cuminc)^2, na.rm = TRUE))
), by = c(
  "time",
  "method",
  "prob_space",
  "failure_time_model",
  "cuminc",
  "censoring_type",
  "X",
  "Z"
)] |>
  ggplot(aes(time, rmse, col = method)) +
  geom_line(aes(linetype = method)) +
  geom_point(aes(shape = method), size = 1.5) +
  facet_grid(
    failure_time_model * censoring_type ~ prob_space,
    labeller = all_labels
  )

preds_df[method != "full"][, .(
  .N,
  rmse = sqrt(mean((pooled_pred - cuminc)^2, na.rm = TRUE))
), by = c(
  "time",
  "method",
  "prob_space",
  "failure_time_model",
  "cuminc",
  "censoring_type",
  "X",
  "Z"
)] |>
  ggplot(aes(time, rmse, col = method)) +
  geom_line(aes(linetype = method)) +
  geom_point(aes(shape = method), size = 1.5) +
  facet_grid(
    failure_time_model * censoring_type ~ prob_space,
    labeller = all_labels
  )





tar_load(true_cuminc_all)

true_cuminc_all[cause == 1][X >= 0 & Z >= 0] |>
  ggplot(aes(time, cuminc)) +
  geom_line(
    aes(
      col = factor(Z),
      linetype = factor(Z)
    )
  ) +
  facet_grid(failure_time_model ~ prob_space * X)
