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

set.seed(558888)
params <- tar_read(true_params_correct_FG_0.15) # but will do for both 0.15 and 0.65
#params <- tar_read(params_weibull_lfps_0.15) # but will do for both 0.15 and 0.65
#tar_load(weibull_FG_lfps_0.15)
# Reduce number of imps n iters for tomorrow..

plan(multisession, workers = 3)
sims <- future_replicate(
  n = 100,
  expr = {
    one_replication(
      args_event_times = list(
        mechanism = "correct_FG",
        #mechanism = "misspec_FG",
        censoring_type = "exponential",
        censoring_params = list("exponential" = "0.49"),
        params = params
      ),
      args_missingness = list(
        mech_params = list(
          "prob_missing" = 0.4,
          #"mechanism_expr" = "1.5 * as.numeric(D == 2)"
          #"mechanism_expr" = "1.5 * scale(log(time))"
          #"mechanism_expr" = "1.5 * ((D == 1) * time + (D != 1) * cens_time)" #now depends on V
          "mechanism_expr" = "1.5 * scale(log(cens_time))"
          # and also one with cens time!!
        )
      ),
      args_imputations = list(
        m = 10, #num_imputations,
        iters = 15, #num_cycles,
        rjlimit = 1000, # not actually needed since X is binary
        rhs_kmi = "1" # NEW switch between 1 vs Z
      ),
      args_predictions = list(timepoints = c(0.001, 0.25, 0.5, 0.75, seq(1, 5, by = 0.5))), #list(timepoints = pred_timepoints),
      true_betas = params$cause1$betas,
      #true_betas = weibull_FG_lfps_0.15[weibull_FG_lfps_0.15[["censoring_type"]] == "exponential", ][["coefs"]]
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

#sims


sims_df <- rbindlist(sims, idcol = "rep_id")
crit <- qnorm((1 - 0.95) / 2, lower.tail = FALSE)

simos <- rbindlist(
  with(
    sims_df,
    Map(
      f = cbind,
      coefs_summary,
      method = method
    )
  ),
  fill = TRUE
)[method != "full"]

simos[, rel_bias := 100 * (estimate - true) / true]

cols <- c("#CABEE9", "#7C7189", "#FAE093", "#D04E59", "#BC8E7D", "#2F3D70")

simos |>
  ggplot(aes(method, rel_bias)) +
  facet_wrap(
    ~ term
  ) +
  coord_flip(ylim = c(-30, 30)) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "none",
    plot.margin = margin(0, 0, 0, 0)
  ) +
  scale_y_continuous(breaks = as.numeric(outer(seq(0, 30, by = 10), c(1, -1)))) +
  scale_x_discrete(limits = rev) +
  geom_hline(
    aes(yintercept = 0),
    linetype = "solid",
    linewidth = 1,
    col = cols[3],
    alpha = 1
  ) +
  stat_summary(
    geom = "segment",
    fun = function(x) mean(x) - crit * (sd(x) / sqrt(length(x))),
    fun.max = function(x) mean(x) + crit * (sd(x) / sqrt(length(x))),
    aes(xend = method, yend = after_stat(ymax), col = method),
    lineend = "round",
    linewidth = 4,
    alpha = 0.75
  ) +
  stat_summary(
    geom = "point", # potentially crossbar
    fun = function(x) mean(x)
    #width = 0.2
  ) +
  scale_color_manual(values = cols[c(1, 2, 6, 4, 5)]) +
  labs(x = "Method", y = "Relative bias in % (Monte Carlo error 95% interval)") +
  theme_light(base_size = 16, base_family = "Roboto Condensed") +
  theme(
    strip.background = element_rect(fill = cols[2], colour = "white"),
    strip.text = element_text(colour = 'white')
  )
