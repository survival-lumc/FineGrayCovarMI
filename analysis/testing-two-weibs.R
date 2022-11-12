library(data.table)
library(ggplot2)
theme_set(theme_bw(base_size = 14))

haz_weib <- function(t, shape, rate, type = "haz") {
  switch(
    type,
    "haz" = shape * rate * t^(shape - 1),
    "cumhaz" = rate * t^shape,
    stop("Try again buddy")
  )
}

# Test function
haz_weib(t = 1:3, shape = 0.5, rate = 0.1)
haz_weib(t = 1:3, shape = 0.5, rate = 0.1, type = "cumhaz")
integrate(haz_weib, lower = 0, upper = 3, shape = 0.5, rate = 0.1)

# Shape of both subdist and cumhaz cause 2
params_sd_cause1 <- c("shape" = 0.75, "rate" = 0.5)
params_cs_cause2 <- c("shape" = 0.75, "rate" = 0.4)

# To check against
integ_fun <- function(t, cause1, cause2) {
  haz_cause1 <- haz_weib(t, cause1[["shape"]], cause1[["rate"]], "haz")
  cumhaz_cause1 <- haz_weib(t, cause1[["shape"]], cause1[["rate"]], "cumhaz")
  cumhaz_cause2 <- haz_weib(t, cause2[["shape"]], cause2[["rate"]], "cumhaz")
  haz_cause1 * exp(-cumhaz_cause1 + cumhaz_cause2)
}

horiz <- 15
integrate(
  f = integ_fun,
  lower = 0,
  upper = horiz,
  cause1 = params_sd_cause1,
  cause2 = params_cs_cause2
)

integ_exact <- function(t, cause1, cause2) {

  a <- cause1[["shape"]]

  # Make sure shapes are the same
  if (!identical(a, cause2[["shape"]])) stop("Shapes not equal!")

  # Simplify if rates are the same - just weibull cumhaz
  b <- cause1[["rate"]]
  if (identical(b, cause2[["rate"]])) {
    haz_weib(t, shape = a, rate = b, type = "cumhaz")
  } else {
    rate_diff <- cause2[["rate"]] - b
    (exp(rate_diff * t^a) - 1) * b / rate_diff
  }
}
integ_exact(horiz, params_sd_cause1, params_cs_cause2)


# Plots -------------------------------------------------------------------


# Recall
params_sd_cause1 <- c("shape" = 0.75, "rate" = 0.5)
params_cs_cause2 <- c("shape" = 0.75, "rate" = 0.4) #0.4, try 0.6

# Check hazards
times <- seq(0.01, 10, by = 0.1)
times_long <- seq(0.01, 100, by = 0.1)
data.table(
  "time" = times,
  "haz_sd" = haz_weib(times, params_sd_cause1[["shape"]], params_sd_cause1[["rate"]]),
  "haz_cs1" = {
    integ_fun(times, params_sd_cause1, params_cs_cause2) /
      (1 - integ_exact(times, params_sd_cause1, params_cs_cause2))
  },
  "haz_cs2" = haz_weib(times, params_cs_cause2[["shape"]], params_cs_cause2[["rate"]])
) |>
  melt.data.table(id.vars = "time", variable.name = "haz_type", value.name = "hazard") |>
  ggplot(aes(time, hazard, col = haz_type)) +
  geom_line(size = 1.25) +
  coord_cartesian(ylim = c(0, 2.5))

# Zoom-in on integral and it's components
data.table(
  "time" = times_long,
  "haz_sd" = haz_weib(times_long, params_sd_cause1[["shape"]], params_sd_cause1[["rate"]]),
  "cumhaz_sd" = haz_weib(times_long, params_sd_cause1[["shape"]], params_sd_cause1[["rate"]], "cumhaz"),
  "cumhaz_cs2" = haz_weib(times_long, params_cs_cause2[["shape"]], params_cs_cause2[["rate"]], "cumhaz"),
  "integral" = integ_exact(times_long, params_sd_cause1, params_cs_cause2),
  "fun_for_integral" = integ_fun(times_long, params_sd_cause1, params_cs_cause2)
) |>
  melt.data.table(id.vars = "time", variable.name = "component", value.name = "hazard") |>
  ggplot(aes(time, hazard, col = component, linetype = "component")) +
  geom_line(size = 1.25) #+
  #coord_cartesian(ylim = c(0, 2.5))

# Try smaller params
params_sd_cause1 <- c("shape" = 0.1, "rate" = 0.5)
params_cs_cause2 <- c("shape" = 0.1, "rate" = 0.1)
plot(1:1000, integ_exact(1:1000, params_sd_cause1, params_cs_cause2))

# Weibul cumhaz goes to infinity!!

# Try gompertz instead??
# See https://sundoc.bibliothek.uni-halle.de/habil-online/07/07H056/t3.pdf

plot(
  1:1000,
  1 - exp(-flexsurv::Hgompertz(1:1000, shape = -0.1, rate = 0.01))
)

# Specify for both



# Gompertz exploration ----------------------------------------------------



library(flexsurv)

# Now Gompertz tings
params_sd_cause1 <- c("shape" = -0.1, "rate" = 0.2)
params_cs_cause2 <- c("shape" = -0.1, "rate" = 0.1) #0.1

1 - exp(-flexsurv::Hgompertz(10000, shape = -0.1, rate = 0.2))

integ_fun <- function(t, cause1, cause2) {
  haz_cause1 <- hgompertz(t, cause1[["shape"]], cause1[["rate"]])
  cumhaz_cause1 <- Hgompertz(t, cause1[["shape"]], cause1[["rate"]])
  cumhaz_cause2 <- Hgompertz(t, cause2[["shape"]], cause2[["rate"]])
  haz_cause1 * exp(-cumhaz_cause1 + cumhaz_cause2)
}

horiz <- 15
integrate(
  f = integ_fun,
  lower = 0,
  upper = horiz, #100000,
  cause1 = params_sd_cause1,
  cause2 = params_cs_cause2
)

fun <- function(t) hgompertz(t, shape = -0.1, rate = 0.2)
fun(horiz)
integrate(fun, lower = 0, upper = horiz)
Hgompertz(horiz, shape = -0.1, rate = 0.2)

integ_exact <- function(t, cause1, cause2) {

  a <- cause1[["shape"]]
  b <- cause1[["rate"]]
  c <- (cause2[["rate"]] - b) / a

  b * (exp(c * (-1 + exp(a * t))) - 1) / (a * c)
}
integ_exact(horiz, cause1 = params_sd_cause1, cause2 = params_cs_cause2)


# Check hazards
times <- seq(0.01, 100, by = 0.1)
times_long <- seq(0.01, 100, by = 0.1)
data.table(
  "time" = times,
  "haz_sd" = hgompertz(times, params_sd_cause1[["shape"]], params_sd_cause1[["rate"]]),
  "haz_cs1" = {
    integ_fun(times, params_sd_cause1, params_cs_cause2) /
      (1 - integ_exact(times, params_sd_cause1, params_cs_cause2))
  },
  "haz_cs2" = hgompertz(times, params_cs_cause2[["shape"]], params_cs_cause2[["rate"]])
) |>
  melt.data.table(id.vars = "time", variable.name = "haz_type", value.name = "hazard") |>
  ggplot(aes(time, hazard, col = haz_type, linetype = haz_type)) +
  geom_line(size = 1.25) +
  coord_cartesian(ylim = c(0, 0.25))


data.table(
  "time" = times_long,
  "haz_sd" = hgompertz(times_long, params_sd_cause1[["shape"]], params_sd_cause1[["rate"]]),
  "cumhaz_sd" = Hgompertz(times_long, params_sd_cause1[["shape"]], params_sd_cause1[["rate"]]),
  "cumhaz_cs2" = Hgompertz(times_long, params_cs_cause2[["shape"]], params_cs_cause2[["rate"]]),
  "integral" = integ_exact(times_long, params_sd_cause1, params_cs_cause2),
  "fun_for_integral" = integ_fun(times_long, params_sd_cause1, params_cs_cause2)
) |>
  melt.data.table(id.vars = "time", variable.name = "component", value.name = "hazard") |>
  ggplot(aes(time, hazard, col = component, linetype = component)) +
  geom_line(size = 1.25)
