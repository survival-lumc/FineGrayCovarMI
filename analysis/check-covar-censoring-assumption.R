#https://github.com/survival-lumc/FineGrayCovarMI/blob/f930d945134da48e59673f9d9749a9b8c0058143/analysis/LSHTM-sims.R

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


# Note for montecarlo tings: use violin instead?? keep monte carlo bits
