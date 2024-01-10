#tar_meta(fields = c(name, time, bytes, seconds, warnings, error)) |> View()


# Overview main sims ------------------------------------------------------


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
censoring_labels <- c(
  "admin" = "Admin. censoring",
  "exponential" = "Exponential censoring",
  "none" = "No censoring"
)
all_labels <- labeller(
  prob_space = pspace_labels,
  censoring_type = censoring_labels,
  failure_time_model = failure_model_labels
)

# Read-in *the* df
simulations_main <- tar_read(simulations_main)

# Unnest
reps_per_scen <- simulations_main[, .(
  coefs = list(rbindlist(coefs_summary, idcol = "rep_id"))
), by = c("prob_space","failure_time_model", "censoring_type", "method")]

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
coefs_main[, ':=' (
  term = ifelse(grepl(pattern = "^X", term), "X", as.character(term)),
  censoring_type = factor(censoring_type, levels = c("none", "exponential", "admin")),
  method = factor(
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
  )
)]


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


# Bias --------------------------------------------------------------------


coefs_main[term == "X" & !(method %in% c("Full data", "Compl. cases"))] |>
  ggplot(aes(method, 100 * (estimate - true) / true)) +
  geom_jitter(aes(col = method), size = 2.5, width = 0.25, alpha = 0.25, shape = 16) +
  facet_grid(
    prob_space ~ failure_time_model * censoring_type,
    labeller = all_labels
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  scale_y_continuous(
    minor_breaks = NULL,
    breaks = c(0, 10, 25, 50, -10, -25, -50)
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
  labs(x = "Method", y = "Relative bias %")


coefs_main[term == "X" & !(method %in% c("Full data", "Compl. cases"))] |>
  ggplot(aes(method, std.error)) +
  geom_jitter(aes(col = method), size = 2.5, width = 0.25, alpha = 0.25, shape = 16) +
  facet_grid(
    prob_space ~ failure_time_model * censoring_type,
    labeller = all_labels
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  stat_summary(
    fun = mean,
    fun.min = mean,
    fun.max = mean,
    geom = "crossbar",
    linewidth = 0.75,
    col = "darkred"
  ) +
  stat_summary(
    data = subset(
      sim_summ$summ,
      stat == "empse" & term == "X" & !(method %in% c("Full data", "Compl. cases"))
    ),
    aes(y = est),
    fun = mean,
    fun.min = mean,
    fun.max = mean,
    geom = "crossbar",
    linewidth = 0.5,
    col = "blue"
  ) +
  scale_color_manual(values = Manu::get_pal("Hoiho")) +
  labs(x = "Method", y = "Standard error")



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
coefs_bis[, "estimate" := estimate - true]

sim_summ <- rsimsum::multisimsum(
  data = coefs_bis,
  par = "term",
  estvarname = "estimate",
  se = "std.error",
  true = 0,
  methodvar = "method",
  ref = "Full data",
  ci.limits = c("conf.low", "conf.high"),
  x = TRUE,
  by = c("censoring_type", "failure_time_model", "prob_space")
)

sim_summ <- simsum(
  data = coefs_bis[term == "X"],
  #par = "term",
  estvarname = "estimate",
  se = "std.error",
  true = 0,
  methodvar = "method",
  ref = "Full data",
  ci.limits = c("conf.low", "conf.high"),
  x = TRUE,
  by = c("censoring_type", "failure_time_model", "prob_space")
)


summary(sim_summ)
autoplot(sim_summ, type = "lolly", par = "X", stats = "bias", )
autoplot(sim_summ, type = "zip", zoom = 0.3)

summary(sim_summ, stats = "cover")


coefs_main[term == "X"][, .(
  n = .N,
  empse = sd(estimate)
), by = c(
  "method",
  "censoring_type",
  "failure_time_model",
  "prob_space"
)] |>
  ggplot(aes(method, empse)) +
  geom_point() +
  facet_grid(
    prob_space ~ failure_time_model * censoring_type,
    labeller = all_labels
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )


coefs_main[term == "X"][, .(
  n = .N,
  coverage = mean(conf.low < true & true < conf.high),
  mcse_covarage = sqrt(
    mean(conf.low < true & true < conf.high) * (1 - mean(conf.low < true & true < conf.high)) / .N
  )
), by = c(
  "method",
  "censoring_type",
  "failure_time_model",
  "prob_space"
)] |>
  ggplot(aes(method, coverage)) +
  geom_point() +
  facet_grid(
    prob_space ~ failure_time_model * censoring_type,
    labeller = all_labels
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )


summary(sim_summ)
autoplot(sim_summ, type = "zip", par = "X")


autoplot(summary(sim_summ), type = "lolly", stats = "bias", par = "X")

autoplot(summary(sim_summ), par = "X", stats = "becover", type = "lolly") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

#autoplot(sim_summ, par = "X", type = "est")

coefs_main[term == "X" & method %in% c("SMC-FCS cause-spec", "SMC-FCS Fine-Gray")] |>
  dcast(
    formula = rep + true + prob_space + failure_time_model + censoring_type ~ method,
    value.var = "estimate"
  ) |>
  ggplot(aes(`SMC-FCS cause-spec` - true, `SMC-FCS Fine-Gray` - true)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  facet_grid(
    prob_space ~ failure_time_model * censoring_type,
    scales = "fixed"
  ) +
  geom_smooth()

coefs_main[term == "X"] |>
  ggplot(aes(method, estimate - true)) +
  geom_jitter(aes(col = method), size = 2.5, width = 0.25, alpha = 0.25, shape = 16) +
  facet_grid(
    prob_space ~ failure_time_model * censoring_type,
    scales = "free_y"
    #labeller = all_labels
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", linewidth = 1) +
  scale_color_manual(values = Manu::get_pal("Hoiho")) +
  labs(x = "Method", y = "Bias")
