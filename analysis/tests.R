library(renv)
#renv::deactivate()

library(data.table)
library(survival)
library(mice)
library(targets)

source("R/data-generation.R")
source("R/kmi-timefixed.R")

set.seed(4894688)
params <- tar_read(true_params_correct_FG_0.15)
dat <- generate_dataset(
  n = 2500,
  args_event_times = list(
    mechanism = "correct_FG",
    params = params,
    censoring_type = "exponential",
    censoring_params = list("exponential" = "0.1 * exp(1 * (as.numeric(X) - 1))")
  ),
  args_missingness = list(
    mech_params = list("prob_missing" = 0.4, "mechanism_expr" = "1.5 * Z")
  )
)
table(dat$D)
plot(survfit(Surv(time, D == 0) ~ X, data = dat))

# Now we add empty vars
dat[, ':=' (
  H_subdist_cause1 = NA_real_,
  newtimes = NA_real_,
  newevent = as.numeric(D == 1)
)]

# Predmats
predmat <- make.predictorMatrix(dat)
predmat[] <- 0
predmat["newtimes", c("time", "D")] <- 1
predmat["H_subdist_cause1", c("newtimes", "D")] <- 1
predmat["X", c("Z", "newevent", "H_subdist_cause1")] <- 1
predmat

impute_cens_times <- function(time, D, X) {

  id_temp <- seq_along(time) # probs can remove
  kmi_single <- kmi(
    Surv(time, D != 0) ~ X,
    data = data.frame(time, D, X),
    etype = D,
    failcode = 1,
    nimp = 1
  )

  imp_dat <- cbind(kmi_single$original.data, kmi_single$imputed.data[[1]])
  new_time <- numeric(length = length(time))
  new_time[D == 2] <- imp_dat[imp_dat$D == 2, ]$newtimes
  new_time[D != 2] <- imp_dat[imp_dat$D != 2, ]$newtimes
  new_time
  # The vector to impute
  #kmi_single$imputed.data[[1]][kmi_single$original.data$D == 2, ]$newtimes
}

# impute_cens_times(dat$time, as.numeric(dat$D) -1, dat$Z)
#
# # something about the confint
# profvis::profvis({
#   kmi_single <- kmi(
#     Surv(time, D != 0) ~ Z,
#     data = data.frame(time = dat$time, D = dat$D, Z = dat$Z),
#     etype = D,
#     failcode = 1,
#     nimp = 1
#   )
# })


# Same as compute marginal cumhaz fun
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


meths <- make.method(dat)
meths["newtimes"] <- paste(
  "~I(", expression(impute_cens_times(time, D, X)),")"
)
meths["H_subdist_cause1"] <- paste(
  "~I(", expression(nelsaalen(newtimes, as.numeric(D == 1))),")"
)

visits <- c("newtimes", "H_subdist_cause1", "X")

imps <- mice(
  data = data.frame(dat),
  method = meths,
  m = 10,
  maxit = 10,
  predictorMatrix = predmat,
  visitSequence = c("newtimes", "H_subdist_cause1", "X")
)
plot(imps)

lapply(complete(imps, "all"), function(x) {
  coxph(Surv(newtimes, newevent) ~ X + Z, data = x)
}) |>
  pool()


# The usual way
kmi_obj <- kmi(
  Surv(time, D != 0) ~ 1,
  data = data.frame(dat[, c("id", "Z", "X", "time", "D", "X_obs")]),
  etype = D,
  failcode = 1,
  nimp = 10
)

impdats_kmi <- lapply(kmi_obj$imputed.data, function(imp_times) {
  impdat <- cbind(kmi_obj$original.data, imp_times)
  impdat$H_subdist_cause1 <- compute_marginal_cumhaz(
    timevar = impdat$newtimes,
    statusvar = impdat$newevent,
    cause = 1,
    type = "cause_spec"
  )
  impdat$newevent <- as.numeric(as.character(impdat$newevent))
  impdat
}
)

impdats_kmi[[1]] |> ncol()
predmat_subdist <- make.predictorMatrix(impdats_kmi[[1]])
meths_subdist <- make.method(impdats_kmi[[1]])
predmat_subdist[] <- 0
predmat_subdist["X", c("Z", "newevent", "H_subdist_cause1")] <- 1

impos_standard <- lapply(impdats_kmi, function(x) {
  impos <- mice(
    data = data.frame(x),
    method = meths_subdist,
    m = 10,
    maxit = 1,
    predictorMatrix = predmat_subdist
  )
  complete(impos, "all")[[1]]
})


lapply(complete(imps, "all"), function(x) {
  coxph(Surv(newtimes, newevent) ~ X + Z, data = x)
}) |>
  pool()

lapply(impos_standard, function(x) {
  coxph(Surv(newtimes, newevent) ~ X + Z, data = x)
}) |>
  pool()

dat[, V := time]
dat[D == 2, V := cens_time]
coxph(Surv(V, D == 1) ~ X_obs + Z, data = dat)

kmi_proper <- kmi(
  Surv(time, D != 0) ~ X_obs, # check this later
  data = data.frame(dat[, c("id", "Z", "X_obs", "time", "D")]),
  etype = D,
  failcode = 1,
  nimp = 10
)
lapply(kmi_proper$imputed.data, function(imp_times) {
  impdat <- cbind(kmi_proper$original.data, imp_times)
  coxph(Surv(newtimes, newevent == 1) ~ X_obs + Z, data = impdat)
}) |>
  pool()

kmi_ign <- kmi(
  Surv(time, D != 0) ~ 1, # check this later
  data = data.frame(dat[, c("id", "Z", "X_obs", "time", "D")]),
  etype = D,
  failcode = 1,
  nimp = 10
)
lapply(kmi_ign$imputed.data, function(imp_times) {
  impdat <- cbind(kmi_ign$original.data, imp_times)
  coxph(Surv(newtimes, newevent == 1) ~ X_obs + Z, data = impdat)
}) |>
  pool()


df1 <- impos_standard[[1]]
df2 <- complete(imps, "all")[[1]]
plot(df1$newtimes, df1$H_subdist_cause1, type = "p")
plot(df1[df1$X == 1, ]$newtimes, df1[df1$X == 1, ]$H_subdist_cause1, type = "p", col = "blue")


plot(df2[df2$X == 0, ]$newtimes, df2[df2$X == 0, ]$H_subdist_cause1, type = "p", col = "blue")

library(ggplot2)
rbindlist(complete(imps, "all"), idcol = ".imp")[D != 2] |>
  ggplot(aes(newtimes, H_subdist_cause1, group = .imp, col = X)) +
  #geom_point() +
  geom_step() #+
  #facet_wrap(~ X, scales = "free")

plot(survfit(Surv(time, D == 0) ~ X_obs, data = dat))

complete(imps, "all")[[1]] |> View()

par(mfrow = c(1, 2))
for (i in seq_len(10)) {
  plot(survfit(Surv(time, D == 0) ~ X, data = impos_standard[[i]]))
  plot(survfit(Surv(time, D == 0) ~ X, data = complete(imps, "all")[[i]]))
}



# Now a different strategy ------------------------------------------------

library(smcfcs)

# Single imp KMI marg
kmi_single <- kmi(
  Surv(time, D != 0) ~ 1, # check this later
  data = data.frame(dat[, c("id", "Z", "X", "time", "D", "X_obs")]),
  etype = D,
  failcode = 1,
  nimp = 1
)
dato <- cbind.data.frame(kmi_single$original.data, kmi_single$imputed.data)
meths_smcfcs <- make.method(dato, defaultMethod = c("norm", "logreg", "mlogit", "podds"))
dato$newevent <- as.numeric(as.character(dato$newevent))

smc_cens <- smcfcs(
  originaldata = dato,
  method = meths_smcfcs,
  smtype = "compet",
  smformula = list(
    "Surv(newtimes, newevent == 1) ~ X + Z",
    "Surv(newtimes, newevent == 0) ~ X"
  ),
  m = 10,
  numit = 30#,
  #rjlimit = 3000
)
plot(smc_cens)

lapply(smc_cens$impDatasets, function(impdat) {
  coxph(Surv(newtimes, newevent == 1) ~ X + Z, data = impdat)
}) |>
  pool()



plot(smc_cens)


smc_ign_cens <- smcfcs(
  originaldata = dato,
  method = meths_smcfcs,
  smtype = "coxph",
  smformula = "Surv(newtimes, newevent) ~ X + Z",
  m = 10,
  numit = 30#,
  #rjlimit = 3000
)
plot(smc_ign_cens)




dat <- generate_dataset(
  n = 2500,
  args_event_times = list(
    mechanism = "correct_FG",
    params = params,
    censoring_type = "exponential",
    censoring_params = list("exponential" = "0.1 * exp(2.5 * (as.numeric(X) - 1))")
  ),
  args_missingness = list(
    mech_params = list("prob_missing" = 0.4, "mechanism_expr" = "1.5 * Z")
  )
)






# Try with admin censoring ------------------------------------------------



dat[, V := time]
dat[D == 2, V := cens_time]
dat[, newevent := as.numeric(D == 1)]
meths_smcfcs <- make.method(dat, defaultMethod = c("norm", "logreg", "mlogit", "podds"))

smc_ign_cens <- smcfcs(
  originaldata = dat,
  method = meths_smcfcs,
  smtype = "coxph",
  smformula = "Surv(V, newevent) ~ X + Z",
  m = 10,
  numit = 30
)
plot(smc_ign_cens)

lapply(smc_ign_cens$impDatasets, function(impdat) {
  coxph(Surv(V, newevent) ~ X + Z, data = impdat)
}) |>
  pool()

smc_cens <- smcfcs(
  originaldata = dat,
  method = meths_smcfcs,
  smtype = "compet",
  smformula = list(
    "Surv(V, D == 1) ~ X + Z",
    "Surv(V, D == 0) ~ X"
  ),
  m = 10,
  numit = 30#,
  #rjlimit = 3000
)
plot(smc_cens)

lapply(smc_cens$impDatasets, function(impdat) {
  coxph(Surv(V, newevent) ~ X + Z, data = impdat)
}) |>
  pool()




# Try simple surv example -------------------------------------------------

library(broom)
library(future)
library(future.apply)

set.seed(9484)

sim_rep <- function() {

  # Gen data
  n <- 2000
  Z <- rnorm(n = n)
  X <- rbinom(n = n, size = 1, prob = 0.3)
  missind_X <- rbinom(n = n, size = 1, prob = 0.4)
  t_tilde <- rexp(n = n, rate = 0.1 * exp(X + Z))
  cens <- rexp(n = n, rate = 0.15 * exp(1 * X))
  time <- pmin(t_tilde, cens)
  D <- as.numeric(t_tilde < cens)
  X_miss <- ifelse(missind_X == 1, NA_real_, X)
  dat <- data.frame(time, D, X = factor(X_miss), X_obs = X, Z)

  # Impute
  meths_smcfcs <- make.method(dat, defaultMethod = c("norm", "logreg", "mlogit", "podds"))

  smc_ign_cens <- smcfcs(
    originaldata = dat,
    method = meths_smcfcs,
    smtype = "coxph",
    smformula = "Surv(time, D) ~ X + Z",
    m = 10,
    numit = 15
  )

  smc_cens <- smcfcs(
    originaldata = dat,
    method = meths_smcfcs,
    smtype = "compet",
    smformula = list(
      "Surv(time, D == 1) ~ X",
      "Surv(time, D == 0) ~ X"
    ),
    m = 10,
    numit = 15#,
    #rjlimit = 3000
  )


  CCA <- coxph(Surv(time, D) ~ X + Z, data = dat)
  CCA_tidy <- tidy(CCA) |> cbind("meth" = "CCA")

  SMC_ign_tidy <- lapply(smc_ign_cens$impDatasets, function(impdat) {
    coxph(Surv(time, D) ~ X + Z, data = impdat)
  }) |>
    pool() |>
    tidy() |>
    cbind(meth = "SMC ignore cens")

  SMC_comp_tidy <- lapply(smc_cens$impDatasets, function(impdat) {
    coxph(Surv(time, D) ~ X + Z, data = impdat)
  }) |>
    pool() |>
    tidy() |>
    cbind(meth = "SMC comp cens")

  rbindlist(list(CCA_tidy, SMC_ign_tidy, SMC_comp_tidy), fill = TRUE)
}

plan(multisession, workers = 3)
test <- future_replicate(
  n = 100,
  expr = sim_rep(),
  simplify = FALSE
)
plan(sequential)

rbindlist(test, idcol = "sim_rep") |>
  ggplot(aes(meth, estimate, fill = meth)) +
  geom_boxplot() +
  geom_jitter(width = 0.15) +
  facet_wrap(~ term) +
  geom_hline(yintercept = 1, linetype = "dotted")+
  stat_summary(fun.y=mean, geom="point", shape=20, size=14, color="red", fill="red")


# Rathouz assumptions



# Hello -------------------------------------------------------------------



set.seed(89465)
n <- 5000
Z <- rnorm(n = n)
X <- rbinom(n = n, size = 1, prob = 0.3)
missind_X <- rbinom(n = n, size = 1, prob = 0.4)
t_tilde <- rexp(n = n, rate = 0.1 * exp(X + Z))
#cens <- rexp(n = n, rate = 0.15 * exp(1 * X))
cens <- rexp(n = n, rate = 0.15 * exp(1 * missind_X))
time <- pmin(t_tilde, cens)
D <- as.numeric(t_tilde < cens)
X_miss <- ifelse(missind_X == 1, NA_real_, X)
dat <- data.frame(time, D, X = factor(X_miss), X_obs = X, Z)

# Impute
meths_smcfcs <- make.method(dat, defaultMethod = c("norm", "logreg", "mlogit", "podds"))

smc_ign_cens <- smcfcs(
  originaldata = dat,
  method = meths_smcfcs,
  smtype = "coxph",
  smformula = "Surv(time, D) ~ X + Z",
  m = 10,
  numit = 15
)

smc_cens <- smcfcs(
  originaldata = dat,
  method = meths_smcfcs,
  smtype = "compet",
  smformula = list(
    "Surv(time, D == 1) ~ X + Z",
    "Surv(time, D == 0) ~ X + Z"
  ),
  m = 10,
  numit = 15#,
  #rjlimit = 3000
)

CCA <- coxph(Surv(time, D) ~ X + Z, data = dat)
CCA_tidy <- tidy(CCA) |> cbind("meth" = "CCA")

SMC_ign_tidy <- lapply(smc_ign_cens$impDatasets, function(impdat) {
  coxph(Surv(time, D) ~ X + Z, data = impdat)
}) |>
  pool() |>
  tidy() |>
  cbind(meth = "SMC ignore cens")

SMC_comp_tidy <- lapply(smc_cens$impDatasets, function(impdat) {
  coxph(Surv(time, D) ~ X + Z, data = impdat)
}) |>
  pool() |>
  tidy() |>
  cbind(meth = "SMC comp cens")

rbindlist(list(CCA_tidy, SMC_ign_tidy, SMC_comp_tidy), fill = TRUE)



# Do another try ----------------------------------------------------------


