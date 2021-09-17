# Testing new function ----------------------------------------------------

# Libraries
library(data.table)
library(prodlim)
library(riskRegression)
library(survival)
library(mice)
library(ggplot2)
library(smcfcs)


# Overall
n <- 2000
predictor_formulas <- list(
  "cause1" = ~ X + Z,
  "cause2" = ~ X + Z
)
X_type = "binary"
source("R/data-generation.R")
source("R/compute-true-cumincs.R")


# Indirect
params_indirect <- list(
  "cause1" = list(
    "betas" = c(0.75, 0.25),
    "base_rate" = 0.15,
    "base_cuminc" = 0.2
  ),
  "cause2" = list(
    "betas" = c(0.5, -0.25),
    "base_rate" = 0.1,
    "base_shape" = 0.5
  )
)

# CS
params_cs <- list(
  "cause1" = list(
    "betas" = c(0.5, 0),
    "base_rate" = 0.1,
    "base_shape" = 0.35
  ),
  "cause2" = list(
    "betas" = c(0.5, -0.25),
    "base_rate" = 0.1,
    "base_shape" = 0.35
  )
)

# Direct
params_direct <- list(
  "cause1" = list(
    "betas" = c(0.5, 0.25),
    "p" = 0.2,
    "base_rate" = 0.5,
    "base_shape" = 1.2
  ),
  "cause2" = list(
    "betas" = c(-0.25, 0.25),
    "base_rate" = 0.15,
    "base_shape" = 0.3
  )
)

# both_fgs
params_fgs <- list(
  "cause1" = list(
    "betas" = c(0.5, 0.25),
    "p" = 0.15,
    "base_rate" = 0.3,
    "base_shape" = 0.6
  ),
  "cause2" = list(
    "betas" = c(-0.5, 0.25),
    "p" = 0.15,
    "base_rate" = 0.3,
    "base_shape" = 0.6
  )
)



# Test data-gen values ----------------------------------------------------

newdat <- data.frame("X" = factor(0, c(0, 1)), "Z" = 0)
t_seq <- seq(0.01, 10, by = 0.1)
true_baseline <- list(
  "indirect" = compute_true_cuminc(t_seq, newdat, params_indirect, "indirect", predictor_formulas),
  "csh_based" = compute_true_cuminc(t_seq, newdat, params_cs, "csh_based", predictor_formulas),
  "direct" = compute_true_cuminc(t_seq, newdat, params_direct, "direct", predictor_formulas),
  "both_fg" = compute_true_cuminc(t_seq, newdat, params_fgs, "both_fg", predictor_formulas)
)

melt(
  rbindlist(true_baseline, idcol = "model_type"),
  id.vars = c("model_type", "time"),
  variable.name = "cause",
  value.name = "cuminc"
) |>
  ggplot(aes(time, cuminc)) +
  geom_line(aes(col = cause, linetype = cause), size = 1.25) +
  facet_wrap(~ model_type, scales = "fixed") +
  coord_cartesian(ylim = c(0, 0.5)) +
  theme_bw()



# Test on actual data -----------------------------------------------------



df_indirect <- generate_complete_dataset(
  n = 2000,
  params = params_indirect,
  model_type = "indirect",
  predictor_formulas,
  X_type = "binary"
)

mod <- FGR(Hist(time, D) ~ X + Z, data = df_indirect, cause = 1)
preds <- predict(mod, newdata = data.frame("X" = factor(0, c(0, 1)), "Z" = 0))
plot(mod$crrFit$uftime, preds, type = "l", ylim = c(0, 0.5))
lines(
  mod$crrFit$uftime,
  compute_true_cuminc(mod$crrFit$uftime, newdat, params_indirect, "indirect", predictor_formulas)$cause1,
  col = "red"
)
