library(data.table)
library(survival)
library(mice)
library(smcfcs)
library(targets)
library(riskRegression)
library(future)
library(future.apply)
library(ggplot2)

source("R/data-generation.R")
source("R/kmi-timefixed.R")
source("R/one-replication.R")

set.seed(4894688)
params <- tar_read(true_params_correct_FG_0.15) # but will do for both 0.15 and 0.65


# Extra sims #1 -----------------------------------------------------------

# - X is MAR-T (depends partially on ID=2), Z is MCAR (to show chained eqs.)
# - correctly specified, random indep cens

dat <- generate_dataset(
  n = 2000,
  args_event_times = list(
    mechanism = "correct_FG",
    params = params,
    censoring_type = "exponential",
    censoring_params = list("exponential" = "0.49 * exp(Z)")
  ),
  args_missingness = list(
    mech_params = list("prob_missing" = 0.4, "mechanism_expr" = "1.5 * time") # important
  )
)
prop.table(table(dat$D))

# Add the mcar bit here, einfach.. but need to change method = FULL


# Extra sims #2 -----------------------------------------------------------


# - p = 0.15 and p = 0.65,
# - correctly specified, X MAR on Z
# - NEW: censoring depending on Z
# - Compare two methods: FG-SMC with and without using Z properly
# .. include other meths?

# One whole rep
plan(multisession, workers = 3)
sims <- future_replicate(
  n = 25,
  expr = {
    one_replication(
      args_event_times = list(
        mechanism = "correct_FG",
        censoring_type = "exponential",
        censoring_params = list("exponential" = "exp(1 * Z)"),
        params = params
        #"correct_FG" = true_params_correct_FG,
        #"misspec_FG" = params_weibull_lfps
      ),
      args_missingness = list(mech_params = list("prob_missing" = 0.4, "mechanism_expr" = "1.5 * Z")),
      args_imputations = list(
        m = 30, #num_imputations,
        iters = 20, #num_cycles,
        rjlimit = 1000, # not actually needed since X is binary
        rhs_kmi = "1" # NEW switch between 1 vs Z
      ),
      args_predictions = list(timepoints = c(0.001, 0.25, 0.5, 0.75, seq(1, 5, by = 0.5))), #list(timepoints = pred_timepoints),
      true_betas = params$cause1$betas,
      #switch(
      #failure_time_model,
      #"correct_FG" = true_params_correct_FG[["cause1"]][["betas"]],
      #"misspec_FG" = weibull_FG_lfps[weibull_FG_lfps[["censoring_type"]] == censoring_type, ][["coefs"]]
      # )
    ) |>
      cbind(prob_space = 0.15)
  },
  simplify = FALSE
)
plan(sequential)

sims_df <- rbindlist(sims, idcol = "rep_id")
crit <- qnorm((1 - 0.95) / 2, lower.tail = FALSE)

rbindlist(
  with(
    sims_df,
    Map(
      f = cbind,
      coefs_summary,
      method = method
    )
  ),
  fill = TRUE
)[method != "full"] |>
  ggplot(aes(method, estimate, fill = method)) +
  geom_violin(trim = T, col = NA) +
  facet_wrap(~ term, scales = "fixed") +
  geom_hline(aes(yintercept = true), linetype = "dotted") +
  # For the monte carlo boxes
  stat_summary(
    fun = function(x) mean(x),
    fun.min = function(x) mean(x) - crit * (sd(x) / sqrt(length(x))),
    fun.max = function(x) mean(x) + crit * (sd(x) / sqrt(length(x))),
    geom = "crossbar",
    alpha = 0.5,
    aes(fill = method),
    fatten = 1,
    linewidth = 0.25
  ) +
  theme_bw() +
  #scale_y_continuous(breaks = seq(0.8, 1.2, by = 0.05)) +
  #coord_cartesian(ylim = c(0.8, 1.2)) +
  labs(y = "Pooled estimate [Monte Carlo SE of bias]")

# Try the pilll
rbindlist(
  with(
    sims_df,
    Map(
      f = cbind,
      coefs_summary,
      method = method
    )
  ),
  fill = TRUE
)[method != "full"] |>
  ggplot(aes(method, estimate, fill = method)) +
  geom_violin(trim = T, col = NA, alpha = 0.5) +
  facet_wrap(~ term, scales = "fixed") +
  geom_hline(aes(yintercept = true), linetype = "dotted") +

  # For the monte carlo boxes
  stat_summary(
    geom = "segment",
    fun = function(x) mean(x) - crit * (sd(x) / sqrt(length(x))),
    fun.max = function(x) mean(x) + crit * (sd(x) / sqrt(length(x))),
    aes(xend = method, yend = after_stat(ymax), col = method),
    lineend = "round",
    linewidth = 3
  ) +
  stat_summary(
    #geom = "crossbar", or just point
    geom = "point",
    fun = function(x) mean(x)
    #width = 0.2
  ) +
  # stat_summary(
  #   fun = function(x) mean(x),
  #   fun.min = function(x) mean(x) - crit * (sd(x) / sqrt(length(x))),
  #   fun.max = function(x) mean(x) + crit * (sd(x) / sqrt(length(x))),
  #   geom = "pointrange",
  #   alpha = 0.5,
  #   aes(col = method),
  #   #fatten = 1,
  #   linewidth = 3#,
  #   #linend = "round"
  # ) +
  theme_minimal() +
  #scale_y_continuous(breaks = seq(0.8, 1.2, by = 0.05)) +
  #coord_cartesian(ylim = c(0.8, 1.2)) +
  labs(y = "Pooled estimate [Monte Carlo SE of bias]")



one_replication(
  args_event_times = list(
    mechanism = "correct_FG",
    censoring_type = "exponential",
    censoring_params = list("exponential" = "exp(1 * Z)"),
    params = params
    #"correct_FG" = true_params_correct_FG,
    #"misspec_FG" = params_weibull_lfps
  ),
  args_missingness = list(mech_params = list("prob_missing" = 0.4, "mechanism_expr" = "1.5 * Z")),
  args_imputations = list(
    m = 30, #num_imputations,
    iters = 20, #num_cycles,
    rjlimit = 1000, # not actually needed since X is binary
    rhs_kmi = "1" # NEW switch between 1 vs Z
  ),
  args_predictions = list(timepoints = c(0.001, 0.25, 0.5, 0.75, seq(1, 5, by = 0.5))), #list(timepoints = pred_timepoints),
  true_betas = params$cause1$betas,
  #switch(
  #failure_time_model,
  #"correct_FG" = true_params_correct_FG[["cause1"]][["betas"]],
  #"misspec_FG" = weibull_FG_lfps[weibull_FG_lfps[["censoring_type"]] == censoring_type, ][["coefs"]]
  # )
) |>
  cbind(prob_space = 0.15)
