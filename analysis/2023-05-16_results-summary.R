# Intermediate summary of results for co-authors 17/5/2023


# General plotting settings -----------------------------------------------


# Set a general theme
invisible(lapply(list.files(here("R"), full.names = TRUE), source))
library(extrafont)
library(Manu)

theme_set(
  theme_light(base_size = 16, base_family = "Roboto Condensed") +
    theme(
      strip.background = element_rect(fill = Manu::get_pal("Hoiho")[[2]], colour = "white"),
      strip.text = element_text(colour = 'white')
    )
)

# Set a bunch of labellers for the plots
pspace_labels <- c("0.15" = "p = 0.15", "0.65" = "p = 0.65")
failure_model_labels <- c("correct_FG" = "Well specified FG", "misspec_FG" = "Misspecified FG")
censoring_labels <- c(
  "curvy_uniform" = "Admin. censoring",
  "exponential" = "Exponential censoring",
  "none" = "No censoring"
)
all_labels <- labeller(
  prob_space = pspace_labels,
  censoring_type = censoring_labels,
  failure_time_model = failure_model_labels
)


# Coefficients (main sims) ------------------------------------------------


coefs_main <- rbindlist(
  with(
    tar_read(simulations_main),
    Map(
      cbind,
      method = method,
      coefs_summary,
      prob_space = prob_space,
      failure_time_model = failure_time_model,
      censoring_type = censoring_type
    )
  )
)

# Let's already read-in the censoring ones
coefs_cens <- rbindlist(
  with(
    tar_read(simulations_cens),
    Map(
      cbind,
      method = method,
      coefs_summary,
      prob_space = prob_space,
      failure_time_model = failure_time_model,
      censoring_type = "exponential",
      cens_rate = cens_rate
    )
  )
)

coefs_all <- rbind(coefs_cens, coefs_main, fill = TRUE)

# Do some cleaning
coefs_all[, term := ifelse(grepl(pattern = "^X", term), "X", as.character(term))]
coefs_all[, method := factor(
  method,
  levels = c("full", "CCA", "mice_comp", "mice_subdist", "smcfcs_comp", "smcfcs_finegray"),
  labels = c(
    "Full data",
    "Compl. cases",
    "MICE cause-spec",
    "MICE subdist",
    "SMC-FCS cause-spec",
    "SMC-FCS Fine-Gray"
  )
)]
coefs_all[, censoring_type := factor(
  censoring_type,
  levels = c("none", "exponential", "curvy_uniform")
)]
coefs_all[is.na(cens_rate) & censoring_type == "exponential", cens_rate := 0.5]
coefs_all[censoring_type == "none", cens_rate := 0]


coefs_all[term == "X" & !(cens_rate %in% c(0.25, 1.5))] |>
  ggplot(aes(method, estimate - true)) +
  geom_jitter(aes(col = method), size = 2.5, width = 0.25, alpha = 0.25, shape = 16) +
  facet_grid(
    prob_space ~ failure_time_model * censoring_type,
    scales = "free_y",
    labeller = all_labels
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", linewidth = 1) +
  stat_summary(
    fun = mean,
    fun.min = mean,
    fun.max = mean,
    geom = "crossbar",
    linewidth = 0.75,
    col = "darkred"
  ) +
  scale_color_manual(values = Manu::get_pal("Hoiho")) +
  labs(x = "Method", y = "Bias")


# Coverage
sim_summ_X <- simsum(
  data = coefs_all[term == "X" & !(cens_rate %in% c(0.25, 1.5))],
  estvarname = "estimate",
  se = "std.error",
  true = "true",
  methodvar = "method",
  ref = "Full data",
  by = c("censoring_type", "failure_time_model", "prob_space")
)

crit <- qnorm((1 - 0.95) / 2, lower.tail = FALSE)

# Coverage
data.table(sim_summ_X$summ)[stat == "cover"] |>
  ggplot(aes(method, est, col = method)) +
  geom_linerange(aes(xmin = method, xmax = method, ymin = 0.6, ymax = est)) +
  geom_point(size = 4) +
  geom_point(
    aes(y = est + crit * mcse),
    shape = 41,
    size = 3
  ) +
  geom_point(
    aes(y = est - crit * mcse),
    shape = 40,
    size = 3
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


coefs_all[term == "X" & !(cens_rate %in% c(0.25, 1.5)), .(
  rmse = sqrt(mean((estimate - true)^2)),
  rmse_mcse = rmse_mcse(estimate, true, .N)
), by = c("failure_time_model", "method", "prob_space", "censoring_type")] |>
  ggplot(aes(method, rmse, col = method)) +
  geom_linerange(aes(xmin = method, xmax = method, ymin = 0, ymax = rmse)) +
  geom_point(size = 4) +
  geom_point(
    aes(y = rmse + crit * rmse_mcse),
    shape = 41,
    size = 3
  ) +
  geom_point(
    aes(y = rmse - crit * rmse_mcse),
    shape = 40,
    size = 3
  ) +
  scale_x_discrete(limits = rev) +
  facet_grid(
    prob_space ~ failure_time_model * censoring_type,
    labeller = all_labels,
    scales = "free_y"
  )  +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    panel.grid.major.y = element_blank()
  ) +
  scale_color_manual(values = Manu::get_pal("Hoiho")) +
  coord_flip(ylim = c(0, 0.25)) +
  labs(y = "RMSE (95% CI with MCSE)", x = "Method")


# Model standard errors? This should look better in a bit..
coefs_all[term == "X" & !(cens_rate %in% c(0.25, 1.5))] |>
  ggplot(aes(method, std.error)) +
  geom_jitter(aes(col = method), size = 2.5, width = 0.25, alpha = 0.25, shape = 16) +
  facet_grid(
    prob_space ~ failure_time_model * censoring_type,
    scales = "free_y",
    labeller = all_labels
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
 # geom_hline(aes(yintercept = 0), linetype = "dashed", size = 1) +
  stat_summary(
    fun = mean,
    fun.min = mean,
    fun.max = mean,
    geom = "crossbar",
    linewidth = 0.75,
    col = "darkred"
  ) +
  scale_color_manual(values = Manu::get_pal("Hoiho")) +
  labs(y = "Model-based standard error X\n(pooled in MI cases)", y = "Bias")



# Investigating censoring -------------------------------------------------


df_cens_plot <- copy(coefs_all[!is.na(cens_rate)])[
  (failure_time_model == "correct_FG" & prob_space == 0.15) |
    (failure_time_model == "misspec_FG" & prob_space == 0.65)
]
df_cens_plot[, cens_rate := factor(
  cens_rate,
  levels = c(0, 0.25, 0.5, 1.5),
  labels = c(0, 0.15, 0.3, 0.5)
)]

df_cens_plot[, .(bias = mean(estimate - true), .N), by = c(
  "method",
  "cens_rate",
  "failure_time_model",
  "term"
)] |>
  ggplot(aes(cens_rate, bias, group = method, col = method)) +
  geom_point(size = 3) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", size = 1) +
  geom_line(aes(linetype = method), size = 1.5) +
  facet_grid(failure_time_model ~ term) +
  scale_color_manual(values = Manu::get_pal("Hoiho")) +
  labs(x = "Proportion censored", y = "Bias")


# Results Z ---------------------------------------------------------------


coefs_main[term == "Z"] |>
  ggplot(aes(method, estimate - true)) +
  geom_jitter(aes(col = method), size = 2.5, width = 0.25, alpha = 0.25, shape = 16) +
  facet_grid(
    prob_space ~ failure_time_model * censoring_type,
    scales = "free_y",
    labeller = all_labels
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", size = 1) +
  stat_summary(
    fun = mean,
    fun.min = mean,
    fun.max = mean,
    geom = "crossbar",
    linewidth = 0.75,
    col = "darkred"
  ) +
  scale_color_manual(values = Manu::get_pal("Hoiho")) +
  labs(x = "Method", y = "Bias")



# Absolute risk preds (main sims) -----------------------------------------


preds_main <- merge(
  tar_read(pooled_preds_main),
  tar_read(true_cuminc_all)
)

# Do some cleaning again
preds_main[, method := factor(
  method,
  levels = c("full", "CCA", "mice_comp", "mice_subdist", "smcfcs_comp", "smcfcs_finegray"),
  labels = c(
    "Full data",
    "Compl. cases",
    "MICE cause-spec",
    "MICE subdist",
    "SMC-FCS cause-spec",
    "SMC-FCS Fine-Gray"
  )
)]
preds_main[, censoring_type := factor(
  censoring_type,
  levels = c("none", "exponential", "curvy_uniform")
)]


# Baseline
preds_main[, .(
  pred = mean(pooled_pred)
), by = c(
  "method",
  "time",
  "prob_space",
  "failure_time_model",
  "censoring_type",
  "cuminc",
  "X",
  "Z"
)][X == 0 & Z == 0] |>
  ggplot(aes(time, 100 * (pred - cuminc))) +
  geom_line(aes(col = method, group = method, linetype = method), size = 2, alpha = 0.7) +
  facet_grid(
    prob_space ~ failure_time_model * censoring_type,
    scales = "free",
    labeller = all_labels
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = Manu::get_pal("Hoiho"))


# X = 1 and Z = 1
preds_main[, .(
  pred = mean(pooled_pred)
), by = c(
  "method",
  "time",
  "prob_space",
  "failure_time_model",
  "censoring_type",
  "cuminc",
  "X",
  "Z"
)][X == 1 & Z == 1] |>
  ggplot(aes(time, 100 * (pred - cuminc))) +
  geom_line(aes(col = method, group = method, linetype = method), size = 2, alpha = 0.7) +
  facet_grid(
    prob_space ~ failure_time_model * censoring_type,
    scales = "free",
    labeller = all_labels
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = Manu::get_pal("Hoiho"))


preds_main[, .(
  pred = mean(pooled_pred)
), by = c(
  "method",
  "time",
  "prob_space",
  "failure_time_model",
  "censoring_type",
  "cuminc",
  "X",
  "Z"
)][X == 1 & Z == 1] |>
  ggplot(aes(time, 100 * pred)) +
  geom_line(aes(col = method, group = method, linetype = method), size = 2, alpha = 0.7) +
  geom_line(aes(y = 100 * cuminc), col = "black", size = 1.25) +
  facet_grid(
    prob_space ~ failure_time_model * censoring_type,
    scales = "free",
    labeller = all_labels
  ) +
  scale_color_manual(values = Manu::get_pal("Hoiho"))


# Big betas example -------------------------------------------------------




coefs_stress <- rbindlist(
  with(
    tar_read(big_betas),
    Map(
      cbind,
      method = method,
      coefs_summary,
      prob_space = prob_space,
      failure_time_model = "correct_FG",
      censoring_type = censoring_type
    )
  )
)

# Do some cleaning
coefs_stress[, term := ifelse(grepl(pattern = "^X", term), "X", as.character(term))]
coefs_stress[, method := factor(
  method,
  levels = c("full", "CCA", "mice_comp", "mice_subdist", "smcfcs_comp", "smcfcs_finegray"),
  labels = c(
    "Full data",
    "Compl. cases",
    "MICE cause-spec",
    "MICE subdist",
    "SMC-FCS cause-spec",
    "SMC-FCS Fine-Gray"
  )
)]
coefs_stress[, censoring_type := factor(
  censoring_type,
  levels = c("none", "exponential", "curvy_uniform")
)]


# Use pmm instead?
coefs_stress |>
  ggplot(aes(method, estimate)) +
  geom_jitter(aes(col = method), size = 2.5, width = 0.25, alpha = 0.5, shape = 16) +
  facet_grid(
    term ~ censoring_type,
    scales = "fixed",
    labeller = all_labels
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  geom_hline(aes(yintercept = true), linetype = "dashed", size = 1) +
  stat_summary(
    fun = mean,
    fun.min = mean,
    fun.max = mean,
    geom = "crossbar",
    linewidth = 0.75,
    col = "darkred"
  ) +
  scale_color_manual(values = Manu::get_pal("Hoiho")) +
  labs(x = "Method", y = "Bias")

# readRDS("data-raw/sims_cens.rds")
# ggsave("interim-sims.png", dpi = 300, units = "in", width = 9, height = 6)


# Horrible zipper plot ----------------------------------------------------



# Shocking attempt at zipper plot
df_coefs[method == "Full data" &
           term == "X" &
           failure_time_model == "correct_FG" &
           prob_space == 0.15 &
           censoring_type == "none", ] |>
  ggplot(aes(estimate - true,
             #cumsum(
             rank(abs((estimate - true) / std.error))
             #  sum(rank(abs((estimate - true) / std.error)))
  )) +
  geom_linerange(aes(xmin = conf.low - true, xmax = conf.high - true))
