invisible(lapply(list.files(here("R"), full.names = TRUE), source))

args_event_times <- list(
  mechanism = "correct_FG",
  #mechanism = "misspec_FG",
  censoring_type = "exponential",
  censoring_params = list(
    "exponential" = 6, #0.5/0.2/0.05 # and also 1.5 = 50% cens;
    "curvy_uniform" = c(0.5, 5),
    "curvyness" = 0.3
  ),
  #params = tar_read(params_weibull_lfps_0.15)
  params = tar_read(true_params_correct_FG_0.15)
)
args_missingness <- list(mech_params = list("prob_missing" = 0, "mechanism_expr" = "Z"))

dat <- generate_dataset(
  n = 2000,
  args_event_times = args_event_times,
  args_missingness = args_missingness
)

table(dat$D)[2]
prop.table(table(dat$D))

testos <- replicate(
  1000, table(generate_dataset(
    n = 2000,
    args_event_times = args_event_times,
    args_missingness = args_missingness
  )$D)[2]
)
hist(testos)

# Censoring checks --------------------------------------------------------




tar_load(all_simulations)
tar_load(cens_sims)

sims_cens_high <- all_simulations[
  censoring_type %in% c("exponential", "none") &
    failure_time_model == "correct_FG" &
    prob_space == 0.15
]
sims_cens_high[, cens_rate := ifelse(censoring_type == "exponential", 0.5, 0)]
cens_tests <- rbind(sims_cens_high, cens_sims, fill = TRUE)

df_coefs_cens <- rbindlist(
  with(
    cens_tests,
    Map(
      cbind,
      method = method,
      coefs_summary,
      cens_rate = cens_rate
    )
  )
)

df_coefs_cens[, term := ifelse(grepl(pattern = "^X", term), "X", as.character(term))]

df_coefs_cens[term == "X"] |>
  #ggplot(aes(method, estimate - true)) +
  ggplot(aes(method, 100 * (estimate - true) / true)) +
  geom_jitter(aes(col = method), size = 2.5, width = 0.25, alpha = 0.5, shape = 16) +
  facet_grid(. ~ cens_rate, scales = "free") +
  theme_bw(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  geom_hline(aes(yintercept = 0), #true),
             linetype = "dashed", size = 1) +
  stat_summary(
    fun = mean,
    fun.min = mean,
    fun.max = mean,
    geom = "crossbar",
    width = 0.75,
    col = "darkred"
  ) +
  scale_color_manual(values = Manu::get_pal("Hoiho")) +
  labs(y = "Relative bias (%)")
  #coord_cartesian(ylim = c(-25, 25))


df_coefs_cens[, .(
  bias = mean(100 * (estimate - true) / true)
), by = c("method", "cens_rate")][, cens_rate := factor(
  cens_rate, levels = c(0, 0.05, 0.2, 0.5, 1.5, 6),
  labels = c(0, 0.05, 0.15, 0.3, 0.5, 0.75)
)] |>
  ggplot(aes(cens_rate, bias, col = method)) +
  geom_hline(yintercept = 0, col = "black", size = 2) +
  geom_line(aes(group = method, linetype = method), size = 2) +
  geom_point(size = 4) +
  theme_bw(base_size = 14) +
  scale_color_manual(values = Manu::get_pal("Hoiho")) +
  labs(y = "Relative bias (%)", x = "Proportion censored") +
  coord_cartesian(ylim = c(-10, 5))



# Test two tings ----------------------------------------------------------




# -- Part 1
args_event_times <- list(
  mechanism = "correct_FG",
  censoring_type = "none",
  params = list(
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
)

args_event_times2 <- list(
  mechanism = "correct_FG",
  censoring_type = "exponential",
  params = list(
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
  ),
  censoring_params = list(
    "exponential" = 1.5, #1.5, #0.2/0.05
    "curvy_uniform" = c(0.5, 5),
    "curvyness" = 0.3
  )
)



args_missingness <- list(mech_params = list("prob_missing" = 0, "mechanism_expr" = "Z"))
args_imputations <- list(m = 10, iters = 20, rjlimit = 1000)
dat <- generate_dataset(
  n = 100000,
  args_event_times = args_event_times,
  args_missingness = args_missingness
)

cs_mod1 <- coxph(Surv(time, D == 1) ~ X + Z, data = dat)
cs_mod2 <- coxph(Surv(time, D == 2) ~ X + Z, data = dat)
coef(cs_mod1)
coef(cs_mod2)
par(mfrow = c(2, 2))
plot(
  cox.zph(cs_mod1, terms = TRUE, transform = "identity"),
  col = "blue", lwd = 3, resid = F,
  ylim = c(-1, 3),
  xlim = c(0, 5)
)
plot(
  cox.zph(cs_mod2, terms = TRUE, transform = "identity"),
  col = "blue", lwd = 3, resid = F,
  ylim = c(-1, 3),
  xlim = c(0, 5)
)



# Experiment --------------------------------------------------------------



args_missingness <- list(mech_params = list("prob_missing" = 0.4,
                                            "mechanism_expr" = "Z"))

dat <- generate_dataset(
  n = 2000,
  args_event_times = args_event_times,
  args_missingness = args_missingness
)

meths_smcfcs <- make.method(dat, defaultMethod = c("norm", "logreg", "mlogit", "podds"))

smcfcs_comp <- smcfcs(
  originaldata = data.frame(dat),
  smtype = "compet",
  smformula = list(
    "Surv(time, D == 1) ~ X + Z",
    "Surv(time, D == 2) ~ X + Z"
  ),
  method = meths_smcfcs,
  rjlimit = args_imputations$rjlimit,
  numit = 20, #args_imputations$iters,
  m = args_imputations$m
)

imp_vals <- rbindlist(smcfcs_comp$impDatasets, idcol = ".imp")[X_missind == 1]

dat_copy <- copy(dat)
dat_copy[, X := X_obs]

rbind(imp_vals, dat_copy, idcol = "imp_or_comp", fill = TRUE)[X_missind == 1] |>
  ggplot(aes(log(time), Z)) +
  geom_point() +
  facet_grid(imp_or_comp ~ D * X) +
  theme_bw(base_size = 14)

dat_copy |>
  ggplot(aes(log(time), Z)) +
  geom_point(data = imp_vals, col = "blue", size = 2, alpha = 0.5,
             shape = "triangle") +
  geom_point(col = "black", alpha = 0.75, size = 2) +

  facet_grid(X_missind ~ D * X) +
  theme_bw(base_size = 14)

