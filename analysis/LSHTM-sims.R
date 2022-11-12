source("R/data-generation.R")

params_direct <- list(
  "cause1" = list(
    "betas" = c(0.5, 0.25),
    "p" = 0.4,
    "base_rate" = 0.5,
    "base_shape" = 1
  ),
  "cause2" = list(
    "betas" = c(0.5, 0.25),
    "base_rate" = 0.15,
    "base_shape" = 1
  )
)

dat <- generate_complete_dataset(
  n = 2000,
  params = params_direct,
  model_type = "direct",
  X_type = "binary", # For now since faster than
  predictor_formulas = list("cause1" = ~  X + Z, "cause2" = ~ X + Z),
  control = list("admin_cens_time" = 1e20, "cens_rate" = 1e-10)
)

dat[]
table(dat$D)

FGR(Hist(time, D) ~ X + Z, data = dat, cause = 1)
coxph(Surv(time, D == 2) ~ X + Z, data = dat)


# Brief sims

map_dfr(
  .x = seq_len(250L),
  .f = ~ {
    dat <- generate_complete_dataset(
      n = 1000,
      params = params_direct,
      model_type = "direct",
      X_type = "continuous", # binary
      predictor_formulas = list("cause1" = ~ X + Z, "cause2" = ~ X + Z),
      control = list("admin_cens_time" = 1e20, "cens_rate" = 1e-10) # add option of NO censoring
    )
    mod_fg <- FGR(Hist(time, D) ~ X + Z, data = dat, cause = 1)$crrFit
    tidy(mod_fg)
    #mod_squeeze <- coxph(Surv(time, D == 2) ~ X + Z, data = dat)
    #bind_rows(tidy(mod_fg), tidy(mod_squeeze), .id = "mod_cause")
  },
  .id = "rep"
) |>
  #mutate(mod_cause = factor(mod_cause, labels = c("FG", "squeeze"))) |>
  ggplot(aes(term, estimate)) +
  geom_boxplot() +
  geom_hline(yintercept = c(0.5, 0.25), linetype = "dotted") #+
  #facet_wrap(~ mod_cause)


#
dat <- generate_complete_dataset(
  n = 50000,
  params = params_direct,
  model_type = "direct",
  X_type = "continuous", # binary
  predictor_formulas = list("cause1" = ~ X + Z, "cause2" = ~ X + Z),
  control = list("admin_cens_time" = 1e20, "cens_rate" = 1e-5) # add option of NO censoring
)
mod_fg <- FGR(Hist(time, D) ~ X + Z, data = dat, cause = 1)$crrFit
tidy(mod_fg)

#other fast packages
library(fastcmprsk)
mod_fg_fast <- fastCrr(Crisk(time, D, cencode = 0) ~ X + Z, variance = TRUE, data = dat)
mod_fg_fast$coef



# Simple MI no censoring --------------------------------------------------

set.seed(984984)
# Perhaps msprep quicker? but FGR much better for prediction

dat <- generate_complete_dataset(
  n = 2000,
  params = params_direct,
  model_type = "direct",
  X_type = "continuous", # binary
  predictor_formulas = list("cause1" = ~ X + Z, "cause2" = ~ X + Z),
  control = list("admin_cens_time" = 1e20, "cens_rate" = 1e-10) # add option of NO censoring
)
table(dat$D)

FGR(Hist(time, D) ~ X + Z, data = dat, cause = 1)
dat[, ':=' (
  D_star = as.numeric(D == 1),
  time_star = ifelse(D == 1, time, max(time) + 1e-8), # works better than setting to super large number..
  miss_ind = rbinom(.N, size = 1, prob = 0.4),
  X_compl = X
)]
coxph(Surv(time_star, D_star) ~ X_compl + Z, data = dat)
dat[, X := ifelse(miss_ind == 1, NA_real_, X)]
make.method(dat, defaultMethod = c("norm", "logreg", "mlogit", "podds"))
imps <- smcfcs(
  originaldata = data.frame(dat),
  smtype = "coxph",
  smformula = "Surv(time_star, D_star) ~ X + Z",
  method = make.method(dat, defaultMethod = c("norm", "logreg", "mlogit", "podds")),
  rjlimit = 5000,
  numit = 30,
  m = 10
)
plot(imps)
lapply(imps$impDatasets, function(imp) coxph(Surv(time_star, D_star) ~ X + Z, data = imp)) |>
  pool() |>
  tidy()

# With censoring ----------------------------------------------------------

set.seed(84645)

dat_cens <- generate_complete_dataset(
  n = 2000,
  params = params_direct,
  model_type = "direct",
  X_type = "continuous", # binary
  predictor_formulas = list("cause1" = ~ X + Z, "cause2" = ~ X + Z),
  control = list("admin_cens_time" = 1e20, "cens_rate" = 0.05) # add option of NO censoring
)

# Test KMI
table(dat_cens$D)
dat_cens[]

kmi_single <- kmi(
  Surv(time, D > 0) ~ 1,
  data = data.frame(dat_cens),
  etype = D,
  failcode = 1,
  nimp = 1
)
# The vector to impute
kmi_single$imputed.data[[1]][kmi_single$original.data$D == 2, ]$newtimes

# Try adjusting
dat_cens[, ':=' (
  t_star = NA_real_,#ifelse(D == 2, NA_real_, time)#,
  cumhaz = NA_real_
  #D_star = as.numeric(D == 1)
)]
dat_cens2 <- copy(dat_cens)
#dat_cens2[, time := NULL]
dat_cens2

meths <- make.method(dat_cens2)
predmat <- make.predictorMatrix(dat_cens2)
predmat[] <- 0
predmat["t_star", c("time", "D")] <- 1
predmat["cumhaz", c("t_star", "D")] <- 1
predmat
#forms <- make.formulas(
#  dat_cens2,

#)
#forms$t_star <- t_star ~ Z

impute_cens_times <- function(time, D) {

  id_temp <- seq_along(time) # probs can remove
  #browser()
  kmi_single <- kmi(
    Surv(time, D > 0) ~ 1,
    data = cbind.data.frame("time" = time, "D" = D),
    etype = D,
    failcode = 1,
    nimp = 1
  )

  imp_dat <- cbind(
    kmi_single$original.data,
    kmi_single$imputed.data[[1]]
  )
  new_time <- numeric(length = length(time))
  new_time[D == 2] <- imp_dat[imp_dat$D == 2, ]$newtimes
  new_time[D != 2] <- imp_dat[imp_dat$D != 2, ]$newtimes
  #new_time[D !=]
  #cbind(
  #  time, D, new_time
  #) |>  View()

  new_time
  # The vector to impute
  #kmi_single$imputed.data[[1]][kmi_single$original.data$D == 2, ]$newtimes
}

nelsaalen <- function(timevar,
                      statusvar,
                      timefix = FALSE) {
  mod <- survival::coxph(
    Surv(time, status) ~ 1,
    control = survival::coxph.control(timefix = timefix),
    data = cbind.data.frame("time" = timevar, "status" = statusvar)
  )
  hazard <- survival::basehaz(mod)
  idx <- match(timevar, hazard[, "time"])
  return(hazard[idx, "hazard"])
}

#nelsaalen(dat_cens3$time, as.numeric(dat_cens3$D == 1))

#impute_cens_times(dat_cens2$time, dat_cens2$D)


meths["t_star"] <- paste(
  "~I(", expression(impute_cens_times(time, D)),")"
)
meths["cumhaz"] <- paste(
  "~I(", expression(nelsaalen(t_star, as.numeric(D == 1))),")"
)
meths
imps <- mice(
  data = data.frame(dat_cens2),
  method = meths,
  m = 10,
  maxit = 25,
  predictorMatrix = predmat
) # next step.. update also nelson aalen

plot(imps)

# Basically no point in updating in terms of nelsaalen; just keep it fixed based on one imp of censoring times
complete(imps, action = "long") |>
  ggplot(aes(t_star, cumhaz, col = factor(.imp), group = factor(.imp))) +
  geom_step() +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(0, 5))# +
  #facet_wrap(~ .imp)

complete(imps, action = "long") |>
  filter(D == 2) |>
  ggplot(aes(t_star, y = factor(.imp), fill = factor(.imp))) +
  ggridges::geom_density_ridges()

imp_dat <- cbind(
  kmi_single$original.data,
  kmi_single$imputed.data[[1]]
)
imp_dat[imp_dat$D == 2, ] |> head()
imp_dat |>  head()

# Generate single KMI

table(dat_cens$D)
# Testing KMI
kmi_obj <- kmi(
  Surv(time, D > 0) ~ 1,
  data = data.frame(dat_cens),
  etype = D,
  failcode = 1,
  nimp = 2
)
kmi_obj

kmi_obj$imputed.data
kmi_obj$imputed.data[[1]]


# Try passive smcfcs. -----------------------------------------------------

dat_cens3 <- copy(dat_cens2)
dat_cens3[, "D_star" := as.numeric(D == 1)]
dat_cens3[, "t_star" := ifelse(D == 2, NA_real_, time)]
dat_cens3[, miss_ind := rbinom(.N, size = 1, prob = 0.4)]
dat_cens3[, X := ifelse(miss_ind == 1, NA_real_, X)]
meths_smcfcs <- make.method(dat_cens3, defaultMethod = c("norm", "logreg", "mlogit", "podds"))
meths_smcfcs
meths_smcfcs["X"] <- "norm"
meths_smcfcs["t_star"] <- "X^2"

smcfcs(
  originaldata = data.frame(dat_cens3),
  smtype = "coxph",
  smformula = "Surv(t_star, D == 1) ~ X + Z",
  method = meths_smcfcs,
  rjlimit = 5000,
  numit = 30,
  m = 10
)

# Will need to again edit source of smcfcs.. ; maybe do for the last meeting!
# Proposal: prepare binary outcome I(D==1) in; but will need to also specify
# comp event indicator somewhere..
# then just prepare new data at the very first step!
# https://github.com/jwb133/smcfcs/blob/00f168de0133b550f5166ba8295f8834d6a960fe/R/smcfcs.r#L454

# To-do: round of KMI, generate all (singly imputed datasets)
# .. then m = 1 in each one

# Capture # of rejection warnings smcfcs
