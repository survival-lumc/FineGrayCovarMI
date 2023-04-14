# Workhorse packages
library("targets")
library("tarchetypes")
library("future")
library("future.callr")
library("here")

# Source other packages/helper functions
source(here("packages.R"))
helper_functions <- list.files(here("R"), full.names = TRUE)
invisible(lapply(helper_functions, source))

# (MAKE INTO TARGETS MARKDOWN!)

# To run pipeline in parallel..
plan(callr)

# For debugging:
#tar_option_set(error = "null")

# Data-generating parameters depend on p-space domination, so we need to iterate
# over that parameter separately
prob_space_domin <- c("low_p" = 0.25, "high_p" = 0.75)

# Rest of scenario parameters are stores in standard data.frame
scenarios <- expand.grid(
  "failure_time_model" = c("correct_FG", "misspec_FG"),
  "censoring_type" = c("none", "exponential"),#, "uniform"),
  stringsAsFactors = FALSE
)

# Here we start the pipeline
simulation_pipeline <- tar_map(
  unlist = FALSE,
  values = list("p" = prob_space_domin),
  # We iterate over the values of p
  tar_target(
    true_params_correct_FG,
    list(
      "cause1" = list(
        "formula" = ~ X + Z,
        "betas" = c(0.75, 0.5),
        "p" = p,
        "base_rate" = 1,
        "base_shape" = 0.75
      ),
      "cause2" = list(
        "formula" = ~ X + Z,
        "betas" = c(0.75, 0.5), # pick some that result in non-zero cause-spec coefs!
        "base_rate" = 1, # might need to increase this to avoid the above?
        "base_shape" = 0.75
      )
    )
  ),
  # Generate large dataset
  tar_target(
    largedat_correct_FG,
    generate_dataset(
      n = 1e3,#1e6,
      args_event_times = list(
        mechanism = "correct_FG",
        params = true_params_correct_FG,
        censoring_type = "none"
      ),
      args_missingness = list(mech_params = list("prob_missing" = 0))
    )
  ),
  # Estimate least-false Weibull parameters, for use in the misspecified FG scenario
  tar_target(
    params_weibull_lfps,
    recover_weibull_lfps(
      large_dat = largedat_correct_FG,
      params_correct_FG = true_params_correct_FG
    )
  ),
  # Generate large dataset using the Weibull parameters, cause-specific hazards
  tar_target(
    largedat_weibull_cause_spec,
    generate_dataset(
      n = 1e3,#1e6,
      args_event_times = list(
        mechanism = "misspec_FG",
        params = params_weibull_lfps,
        censoring_type = "none"
      ),
      args_missingness = list(mech_params = list("prob_missing" = 0))
    )
  ),
  # Now recover the least-false weibull parameters for that setting,
  # .. these are our 'true' values
  tar_target(
    largedat_weibull_FG_lfps,
    recover_fg_lps(largedat_weibull_cause_spec)
  ),

  # This are the actual simulation replication, iterate over the remaining scenario parameters
  tar_map_rep(
    name = prob_space, # for the nice names after
    values = scenarios,
    command = one_replication(
      args_event_times = list(
        mechanism = failure_time_model,
        params = switch(
          failure_time_model,
          "correct_FG" = true_params_correct_FG,
          "misspec_FG" = params_weibull_lfps
        ),
        censoring_type = censoring_type
      ),
      args_missingness = list(
        mech_params = list("prob_missing" = 0.1, "mechanism_expr" = "Z") # remember to change to vals in sim protocol
      ),
      args_imputations = list(m = 2, iters = 1, rjlimit = 1000)
    ),
    reps = 2,
    batches = 1,
    combine = TRUE
  )
)

# Here we bring together all the simulation scenarios
list(
  simulation_pipeline,
  tar_combine(
    all_sims,
    simulation_pipeline[["prob_space"]],
    command = dplyr::bind_rows(!!!.x, .id = "prob_space")
  )
)

# Might need to add least-false true onto the dataframe
