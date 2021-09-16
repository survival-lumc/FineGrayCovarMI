time_grid <- seq(0.1, 10, by = 0.1)



true_cuminc_indirect <- function(t,
                                 params_sd_cause1,
                                 params_cs_cause2,
                                 x) {
  # Cause 1 first
  hr_sd <- exp(params_sd_cause1[["beta_X"]] * x[["X"]] + params_sd_cause1[["beta_Z"]] * x[["Z"]])
  cumhaz_cause1 <- haz_sd_cause1(t = t, x = x, params_sd_cause1, type = "cumulative")
  cuminc_cause1 <- 1 - exp(-cumhaz_cause1 * hr_sd)

  # Cause 2
  prod <-  function(t) {
    haz <- haz_cs_cause2(t, x, params_cs_cause2, type = "hazard")
    cumhaz_cause1 <- haz_cs_cause1(t, x, params_sd_cause1, params_cs_cause2, type = "cumulative")
    cumhaz_cause2 <- haz_cs_cause2(t, x,params_cs_cause2, type = "cumulative")

    haz * exp(-cumhaz_cause1 - cumhaz_cause2)
  }

  ci_func <- Vectorize(function(upp) {
    if (upp == 0) upp <- .Machine$double.eps
    integrate(prod, lower = 0, upper = upp)$value
  })

  res <- cbind.data.frame(
    "time" = t,
    "cause1" = cuminc_cause1,
    "cause2" = ci_func(t)
  )

  return(res)
}

df_cs <- true_cuminc_indirect(time_grid, params_sd_cause1, params_cs_cause2, list("X" = 0, "Z" = 0))

plot(
  x = time_grid,
  y = df_cs$cause1,
  type = "l",
  ylim = c(0, 0.6)
)
lines(time_grid, df_cs$cause2, col = "blue")

true_cuminc_cs <- function(t, params_cs_cause1, params_cs_cause2, x) {

  prod <-  function(t, cause) {

    haz <- if (cause == 1) {
      haz_cs_cause2(t, x, params_cs_cause1, type = "hazard") # same func
    } else {
      haz_cs_cause2(t, x, params_cs_cause2, type = "hazard")
    }

    cumhaz_cause1 <- haz_cs_cause2(t, x, params_cs_cause1, type = "cumulative")
    cumhaz_cause2 <- haz_cs_cause2(t, x, params_cs_cause2, type = "cumulative")

    haz * exp(-cumhaz_cause1 - cumhaz_cause2)
  }

  ci_func <- Vectorize(function(upp, ...) {
    if (upp == 0) upp <- .Machine$double.eps
    integrate(prod, lower = 0, upper = upp, ...)$value
  })

  res <- cbind.data.frame(
    "time" = t,
    "cause1" = ci_func(t, cause = 1),
    "cause2" = ci_func(t, cause = 2)
  )

  return(res)
}

df_cs <- true_cuminc_cs(time_grid, params_cs_cause1, params_cs_cause2, list("X" = 0, "Z" = 0))

plot(
  x = time_grid,
  y = df_cs$cause1,
  type = "l",
  ylim = c(0, 0.6)
)
lines(time_grid, df_cs$cause2, col = "blue")



# Direct stuff ------------------------------------------------------------



true_cuminc_direct <- function(t, params_cuminc_cause1, params_condit_cause2, x) {

  # Cause 1
  hr_sd <- exp(
    params_cuminc_cause1[["beta_X"]] * x[["X"]] + params_cuminc_cause1[["beta_Z"]] * x[["Z"]]
  )
  p <- params_cuminc_cause1[["p"]]
  shape <- params_cuminc_cause1[["base_shape"]]
  rate <- params_cuminc_cause1[["base_rate"]]

  cuminc_cause1 <- 1 - (1 - p * (1 - exp(-rate * t^shape)))^hr_sd

  # Cause 2
  px_2 <- (1 - p)^hr_sd
  hr_condit <- exp(
    params_condit_cause2[["beta_star_X"]] * x[["X"]] + params_condit_cause2[["beta_star_Z"]] * x[["Z"]]
  )
  cumhaz_condit <- params_condit_cause2[["base_rate"]] * hr_condit * t^params_condit_cause2[["base_shape"]]
  cuminc_cause2 <- (1 - exp(-cumhaz_condit)) * px_2

  res <- cbind.data.frame(
    "time" = t,
    "cause1" = cuminc_cause1,
    "cause2" = cuminc_cause2
  )

  return(res)
}


df_direct <- true_cuminc_direct(time_grid, params_cuminc_cause1,
                                params_condit_cause2, list("X" = 0, "Z" = 0))

plot(
  x = time_grid,
  y = df_direct$cause1,
  type = "l",
  ylim = c(0, 0.6)
)
lines(time_grid, df_direct$cause2, col = "blue")


# FGs ---------------------------------------------------------------------


true_cuminc_fgs <- function(t, params_cuminc_cause1, params_cuminc_cause2, x) {

  # Cause 1
  hr_sd1 <- exp(
    params_cuminc_cause1[["beta_X"]] * x[["X"]] + params_cuminc_cause1[["beta_Z"]] * x[["Z"]]
  )
  p1 <- params_cuminc_cause1[["p"]]
  shape1 <- params_cuminc_cause1[["base_shape"]]
  rate1 <- params_cuminc_cause1[["base_rate"]]

  cuminc_cause1 <- 1 - (1 - p1 * (1 - exp(-rate1 * t^shape1)))^hr_sd1

  # Cause 2
  hr_sd2 <- exp(
    params_cuminc_cause2[["beta_X"]] * x[["X"]] + params_cuminc_cause2[["beta_Z"]] * x[["Z"]]
  )
  p2 <- params_cuminc_cause2[["p"]]
  shape2 <- params_cuminc_cause2[["base_shape"]]
  rate2 <- params_cuminc_cause2[["base_rate"]]

  cuminc_cause2 <- 1 - (1 - p2 * (1 - exp(-rate2 * t^shape2)))^hr_sd2

  res <- cbind.data.frame(
    "time" = t,
    "cause1" = cuminc_cause1,
    "cause2" = cuminc_cause2
  )

  return(res)
}

df_fgs <- true_cuminc_fgs(time_grid, params_cuminc_cause1,
                             params_cuminc_cause2, list("X" = 1, "Z" = 3))

plot(
  x = time_grid,
  y = df_fgs$cause1,
  type = "l",
  ylim = c(0, 0.6)
)
lines(time_grid, df_fgs$cause2, col = "blue")
