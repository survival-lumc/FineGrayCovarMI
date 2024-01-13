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
  scale_y_continuous(minor_breaks = NULL, breaks = c(0, 10, 25, 50, 75, -10, -25, -50, -75)) +
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
dodge_w <- 0.75

df_summ <- data.table(sim_summ$summ)
df_summ_var <- df_summ[stat %in% c("modelse", "empse", "cover") & method != "Full data" & term == "X"]
df_summ_var[, stat := factor(
  stat, levels = c("empse", "modelse", "cover"),
  labels = c("Emp. SE", "Mod. SE", "Coverage")
)]
df_summ_var[, stat_lab := paste0(
  "bold('", stat, ":')~",
  round(est, 3), "~(", round(mcse, 2), ")"
)]


df_summ_var[stat %in% c("Emp. SE", "Mod. SE")] |>
  ggplot(aes(method, est, col = method)) +
  geom_text(
    data = df_summ_var[stat == "Coverage"],#merge(df_lab, left_lab),
    aes(label = stat_lab, y = 0.21),
    hjust = 0,
    parse = TRUE,
    family = "Roboto Condensed",
    size = 3
  ) +
  geom_linerange(
    aes(xmin = method, xmax = method, ymin = 0, ymax = est, linetype = stat),
    position = position_dodge(width = dodge_w),
    linewidth = 0.5
  ) +
  geom_point(
    size = 2,
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
  coord_flip(ylim = c(0.05, 0.31)) +
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
  scale_color_manual(values = cols[c(1, 2, 6, 4, 5)]) +
  labs(x = "Method", y = "Standard error (95% Monte-Carlo error interval)")

ggsave(
  filename = "analysis/figures/perf_X.pdf",
  width = 7,
  scale = 1.25,
  height = 10,
  device = cairo_pdf
)



# Trying prediction plots -------------------------------------------------


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

# For now just baseline
preds_main[is.na(imp), imp := 0]
new_pat <- c(0, 0) # Then c(1, 1)

# Calculate linear predictor (do this bit in targets later)
# .. for baseline this bit is unnecessary!
# preds_main[, "lp" := sum(unlist(coefs) * new_pat), by = c(
#   "method",
#   "imp",
#   "time",
#   "rep_id",
#   "prob_space",
#   "failure_time_model",
#   "censoring_type"
# )]
preds_main[, lp := 0]

# Calculate cumulative incidence
preds_main[, pred := 1 - (1 - base_cuminc)^exp(lp)]

# Pool after cloglog!
pooled_preds <- preds_main[, .(pooled_pred = inv_cloglog(mean(cloglog(pred)))), by = c(
  "method",
  "rep_id",
  "time",
  "prob_space",
  "failure_time_model",
  "censoring_type"
)]
pooled_preds[, ':=' (X = 0, Z = 0)]

# Merge them with the true values
true_vals <- tar_read(true_cuminc_all)[cause == 1]
preds_df <- merge(pooled_preds, true_vals)

# Again we can do the measures in the same way
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

#preds_df[is.na(pooled_pred)] # only one rep had NA predictions
df_summ_pred <- data.table(sim_summ_preds$summ)[censoring_type != "admin"]
df_summ_pred[, time := as.numeric(as.character(time))]
df_summ_pred$stat |> unique() # no SEs recorded so just empirical ones

# Let's make some plots - focus only on no and random cens (since admin not much going on)
# And we have to keep method = FULL so confirm it is mainly model misspec

p_cumincs <- df_summ_pred[stat == "thetamean"] |>
  ggplot(aes(time, est, group = method, col = method)) +
  geom_line(aes(linetype = method), linewidth = 1) +
  geom_point(aes(shape = method)) +
  facet_grid(
    prob_space * failure_time_model * censoring_type ~ .,
    labeller = all_labels,
    scales = "free"
  ) +
  geom_line(
    data = preds_df[censoring_type != "admin"],
    aes(y = cuminc),
    col = "black",
    linewidth = 1
  ) +
  scale_color_manual(values = cols[c(3, 1, 2, 6, 4, 5)]) +
  theme(
    #legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_blank()
  ) +
  labs(
    col = "Method",
    linetype = "Method",
    shape = "Method",
    x = "Time",
    y = "Baseline cumulative incidence cause 1"
  )


p_rmse <- df_summ_pred[stat == "mse"] |>
  ggplot(aes(time, sqrt(est), group = method, col = method)) +
  geom_line(aes(linetype = method), linewidth = 1) +
  geom_point(aes(shape = method)) +
  facet_grid(
    prob_space * failure_time_model * censoring_type ~ .,
    labeller = all_labels,
    scales = "fixed" # free?
  ) +
  scale_color_manual(values = cols[c(3, 1, 2, 6, 4, 5)]) +
  labs(
    col = "Method",
    linetype = "Method",
    shape = "Method",
    x = "Time",
    y = "Root mean square error (RMSE)"
  )

comb <- p_cumincs + p_rmse & theme(legend.position = "top")
comb + plot_layout(guides = "collect")

ggsave(
  filename = "analysis/figures/preds_baseline.pdf",
  width = 7,
  scale = 1.5,
  height = 10,
  device = cairo_pdf
)


# Now with X = 1 and  Z =1  -----------------------------------------------


new_pat <- c(1, 1) # Then c(1, 1)

# Calculate linear predictor (do this bit in targets later)
# .. for baseline this bit is unnecessary!
preds_main[, "lp" := sum(unlist(coefs) * new_pat), by = c(
  "method",
  "imp",
  "time",
  "rep_id",
  "prob_space",
  "failure_time_model",
  "censoring_type"
)]

# Calculate cumulative incidence
preds_main[, pred := 1 - (1 - base_cuminc)^exp(lp)]

# Pool after cloglog!
pooled_preds <- preds_main[, .(pooled_pred = inv_cloglog(mean(cloglog(pred)))), by = c(
  "method",
  "rep_id",
  "time",
  "prob_space",
  "failure_time_model",
  "censoring_type"
)]
pooled_preds[, ':=' (X = 1, Z = 1)]

# Merge them with the true values
true_vals <- tar_read(true_cuminc_all)[cause == 1]
preds_df <- merge(pooled_preds, true_vals)

# Again we can do the measures in the same way
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

#preds_df[is.na(pooled_pred)] # only one rep had NA predictions
df_summ_pred <- data.table(sim_summ_preds$summ)[censoring_type != "admin"]
df_summ_pred[, time := as.numeric(as.character(time))]

# Let's make some plots - focus only on no and random cens (since admin not much going on)
# And we have to keep method = FULL so confirm it is mainly model misspec

p_cumincs <- df_summ_pred[stat == "thetamean"] |>
  ggplot(aes(time, est, group = method, col = method)) +
  geom_line(aes(linetype = method), linewidth = 1) +
  geom_point(aes(shape = method)) +
  facet_grid(
    prob_space * failure_time_model * censoring_type ~ .,
    labeller = all_labels,
    scales = "free"
  ) +
  geom_line(
    data = preds_df[censoring_type != "admin"],
    aes(y = cuminc),
    col = "black",
    linewidth = 1
  ) +
  scale_color_manual(values = cols[c(3, 1, 2, 6, 4, 5)]) +
  theme(
    #legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_blank()
  ) +
  labs(
    col = "Method",
    linetype = "Method",
    shape = "Method",
    x = "Time",
    y = "Cumulative incidence cause 1 (X = 1 and Z = 1)"
  )


p_rmse <- df_summ_pred[stat == "mse"] |>
  ggplot(aes(time, sqrt(est), group = method, col = method)) +
  geom_line(aes(linetype = method), linewidth = 1) +
  geom_point(aes(shape = method)) +
  facet_grid(
    prob_space * failure_time_model * censoring_type ~ .,
    labeller = all_labels,
    scales = "fixed" # free?
  ) +
  scale_color_manual(values = cols[c(3, 1, 2, 6, 4, 5)]) +
  labs(
    col = "Method",
    linetype = "Method",
    shape = "Method",
    x = "Time",
    y = "Root mean square error (RMSE)"
  )

comb <- p_cumincs + p_rmse & theme(legend.position = "top")
comb + plot_layout(guides = "collect")

ggsave(
  filename = "analysis/figures/preds_X1Z1.pdf",
  width = 7,
  scale = 1.5,
  height = 10,
  device = cairo_pdf
)



# What are the true ones? -------------------------------------------------


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
