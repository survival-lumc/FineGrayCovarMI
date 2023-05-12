# Tests for number of imputations
library(targets)
library(here)
library(smcfcs)
library(mice)
library(data.table)
library(survival)
library(kmi)
library(ggplot2)
invisible(lapply(list.files(here("R"), full.names = TRUE), source))

args_event_times <- list(
  #mechanism = "correct_FG",
  mechanism = "misspec_FG",
  censoring_type = "exponential",
  #params = tar_read(params_weibull_lfps_0.15)
  params = tar_read(true_params_correct_FG_0.15)
)
args_missingness <- list(mech_params = list("prob_missing" = 0.4, "mechanism_expr" = "Z"))
args_imputations <- list(m = 2, iters = 20, rjlimit = 1000)

dat <- generate_dataset(
  n = 2000,
  args_event_times = args_event_times,
  args_missingness = args_missingness
)
dat[, D := as.numeric(as.character(D))]
dat
meths_smcfcs <- make.method(dat, defaultMethod = c("norm", "logreg", "mlogit", "podds"))

table(dat$D)

imps <- smcfcs.finegray(
  originaldata = data.frame(dat),
  smformula = "Surv(time, D) ~ X + Z",
  method = meths_smcfcs,
  cause = 1,
  m = 20,
  numit = 100,
  rjlimit = 1000
)

plot(imps, nrow = 2) +
  theme_light() +
  theme(legend.position = "none") +
  #geom_line(linewidth = 1.5, alpha = 0.25) +
  scale_x_continuous(breaks = seq(0, 100, by = 10))

imp_mods <- lapply(
  imps$impDatasets,
  function(imp) coxph(Surv(newtimes, newevent) ~ X + Z, data = imp, x = TRUE)
)

howManyImputations::how_many_imputations(imp_mods, )

libr



p <- plot(imps)

library(ggmcmc)
library(coda)
library(tidyverse)
df <- p$data |>
  pivot_wider(id_cols = c(iters, imp), names_from = covar)

obj <- lapply(
  split(df, df$imp),
  function(x) as.mcmc(x[, c("X1", "Z")])
) |>
  as.mcmc.list() |>
  ggs()

ggs_traceplot(obj)
ggs_running(obj)


# Check proportionality ---------------------------------------------------



args_event_times <- list(
  mechanism = "correct_FG",
  #mechanism = "misspec_FG",
  censoring_type = "exponential",
  #params = tar_read(params_weibull_lfps_0.15)
  params = tar_read(true_params_correct_FG_0.15)
)
args_missingness <- list(mech_params = list("prob_missing" = 0, "mechanism_expr" = "Z"))

dat <- generate_dataset(
  n = 100000,
  args_event_times = args_event_times,
  args_missingness = args_missingness
)

coxph(Surv(time, D == 0) ~ X + Z, data = dat)
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


#plot(cox.zph(cs_mod1, terms = TRUE), col = "blue", lwd = 3, resid = F, ylim = c(-0.5, 3))
#plot(cox.zph(cs_mod2, terms = TRUE), col = "blue", lwd = 3, resid = F, ylim = c(-0.5, 3))

# Test now plot no censoring scenario -------------------------------------


invisible(lapply(list.files(here("R"), full.names = TRUE), source))

generate_covariates <- function(n) {
  dat <- data.table(id = seq_len(n), Z = rnorm(n = n, mean = 0, sd = 1))
  dat[, X := rnorm(n, mean = 0)] #  indep
  return(dat)
}


# TRY THIS IN THE MISSPECIFIED CASE
args_event_times <- list(
  mechanism = "correct_FG",
  #mechanism = "misspec_FG",
  censoring_type = "exponential",
  #censoring_type = "none",
  #params = tar_read(params_weibull_lfps_0.15)
  #params = tar_read(true_params_correct_FG_0.15)
  params = list(
    "cause1" = list(
      "formula" = ~ X + Z,
      "betas" = c(0.75, 0),
      "p" = 0.15,
      "base_rate" = 1,
      "base_shape" = 0.75
    ),
    "cause2" = list(
      "formula" = ~ X + Z,
      "betas" = c(0.75, 0),
      "base_rate" = 1,
      "base_shape" = 0.75
    )
  ),
  censoring_params = list(
    "exponential" = 0.05, #0.2/0.05
    "curvy_uniform" = c(0.5, 5),
    "curvyness" = 0.3
  )
)
args_missingness <- list(mech_params = list("prob_missing" = 0.4, "mechanism_expr" = "Z"))
args_imputations <- list(m = 10, iters = 20, rjlimit = 1000)

dat <- generate_dataset(
  n = 2000,
  args_event_times = args_event_times,
  args_missingness = args_missingness
)

#plot(prodlim(Hist(time, D) ~ 1, data = dat), cause = 1)
#plot(prodlim(Hist(time, D) ~ 1, data = dat), cause = 2, add = TRUE)
coxph(Surv(time, D == 1) ~ X + Z, data = dat) |> coef()
coxph(Surv(time, D == 2) ~ X + Z, data = dat) |> coef()

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
  numit = 50, #args_imputations$iters,
  m = args_imputations$m
)

plot(smcfcs_comp) +
  theme_light() +
  theme(legend.position = "none") +
  geom_line(linewidth = 1) +
  scale_x_continuous(breaks = seq(0, 50, by = 5)) +
  coord_cartesian(ylim = c(0, 1.5))

# Investigate imps
imps_long <- dplyr::bind_rows(smcfcs_comp$impDatasets, .id = ".imp")
dat$`.imp` <- 0L
imps_long_complete <- rbind.data.frame(dat, imps_long)
smcfcs_mids <- as.mids(long = imps_long_complete, .id = "id")
mice::stripplot(smcfcs_mids)


mods_smcfcs <- pbapply::pblapply(
  smcfcs_comp$impDatasets,
  function(imp) FGR(Hist(time, D) ~ X + Z, data = imp, cause = 1)$crrFit
)

pool(mods_smcfcs) |> tidy()



# COnditional dist checking
t <- seq(0, 10, by = 0.1)
haz_subdist <- function(t, a, b, p) {
  nom <- p * a * b * exp(-b * t^a) * t^(a - 1)
  denom <- 1 - p * (1 - exp(-b * t^a))
  nom / denom
}

cumhaz_subdist <- Vectorize(function(t, a, b, p) {
  if (t > 0) {
    integrate(
      haz_subdist,
      a = a,
      b = b,
      p = p,
      lower = 0L,
      upper = t
    )$value
  } else 0
})


# It fookin works!!
plot(t, haz_subdist(t, a = 0.75, b = 1, p = 0.65))
plot(t, cumhaz_subdist(t, a = 0.75, b = 1, p = 0.15))
plot(
  t,
  1 - exp(-cumhaz_subdist(t, a = 1.5, b = 1, p = 0.15)),
  ylim = c(0, 0.8)
)


# Let's try adding to data
dat <- generate_dataset(
  n = 100000,
  args_event_times = args_event_times,
  args_missingness = args_missingness
)
dat[, cumhaz_subdist_base_true := cumhaz_subdist(
  t = time,
  a = args_event_times$params$cause1$base_shape,
  b = args_event_times$params$cause1$base_rate,
  p = args_event_times$params$cause1$p
)]

dat |>
  ggplot() + #, linetype = factor(Z))) +
  geom_smooth(
    aes(x = cumhaz_subdist_base_true, y = X_obs, col = D),
    linetype = "dotted",
    method = "gam",
    se = FALSE,
    formula = y ~ s(x, bs = "cs"),
    size = 1.5
  ) +
  scale_color_manual(
    "D",
   values = Manu::get_pal("Hoiho")#,
  #labels = c("Censored", "Observed")
  ) +
  # guides(lty = guide_legend("Z", override.aes = list(col = "black"))) +
  labs(x = "Estimated Cumulative subdist baseline hazard at T") +
  theme_bw(base_size = 14) +
  coord_cartesian(xlim = c(0, 0.15), ylim = c(-1.5, 0.75))


dat[, H1 := compute_marginal_cumhaz(
  timevar = time,
  statusvar = D,
  cause = 1,
  type = "cause_spec"
)]
dat[, H2 := compute_marginal_cumhaz(
  timevar = time,
  statusvar = D,
  cause = 2,
  type = "cause_spec"
)]


dat |>
  ggplot() + #, linetype = factor(Z))) +
  geom_smooth(
    aes(x = H1, y = X_obs, col = D),
    linetype = "dotted",
    method = "gam",
    se = FALSE,
    formula = y ~ s(x, bs = "cs"),
    size = 1.5
  ) +
  geom_smooth(
    aes(x = H2, y = X_obs, col = D),
    linetype = "dashed",
    method = "gam",
    se = FALSE,
    formula = y ~ s(x, bs = "cs"),
    size = 1.5
  ) +
  #scale_color_manual(
  #  "D",
   # values = Manu::get_pal("Hoiho")#,
    #labels = c("Censored", "Observed")
  #) +
  # guides(lty = guide_legend("Z", override.aes = list(col = "black"))) +
  labs(x = "Estimated Cumulative subdist baseline hazard at T") #+
  #coord_cartesian(xlim = c(0, 5))
