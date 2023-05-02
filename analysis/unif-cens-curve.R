dat <- generate_dataset(
  n = 1e3,
  list(
    mechanism = "correct_FG",
    params = tar_read(true_params_correct_FG_0.25),
    censoring_type = "none"#,
    #censoring_params = list("exponential" = 0.45, "uniform" = c(0, 3))
  ),
  list(mech_params = list("prob_missing" = 0, "mechanism_expr" = "1"))
)


# The curvy uniform

# We want admin cens on 6 months (so 0.5 years) to 5 years
# Uniform on the period has density 1 / (5 - 0.5)
x <- seq(0.5, 5, by = 0.01)
plot(x, dunif(x, 0.5, 5), type = "l")

# And corresponding CDF, with equation (x - 0.5) / (5 - 0.5) = (x - 0.5) / 4.5
# (need integration constant)
plot(x, punif(x, 0.5, 5), type = "l")

# What we want is to add a little curvature to this cdf, so that we censor
#.. some more people early on. Easy way: raise the CDF to a fractional power
#.. this dist will still be super easy to sample from

curvyness <- 0.25
plot(x, (x - 0.5) / 4.5, type = "l")
lines(x, ((x - 0.5) / 4.5)^curvyness, lty = "dashed")

# Can also plot the density, which is just the derivative (to-do)

# Sample censoring times
cens_times <- 4.5 * runif(1e4)^(1 / curvyness) + 0.5
hist(cens_times, breaks = 50, freq = FALSE)

# Hurray!
