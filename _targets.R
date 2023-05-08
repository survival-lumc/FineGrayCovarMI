# Possible to-do/ideas:
# - record rej sampling warnings if using contin X?
# - Number of imputations per kmi dataset change?
# - Add proposed smcfcs.finegray to smcfcs() package
# - Visualise missing mechanism with jitter plots! (visualise the settings, together with
# base cumincs)

# Workhorse packages
library("targets")
library("tarchetypes")
library("future")
library("future.callr")
library("here")

# Source other packages/helper functions
source(here("packages.R"))
invisible(lapply(list.files(here("R"), full.names = TRUE), source))

# (MAKE INTO TARGETS MARKDOWN!)

# To run pipeline in parallel..
plan(callr)

# For debugging:
#tar_option_set(error = "null")

# Data-generating parameters depend on p-space domination, so we need to iterate
# over that parameter separately
prob_space_domin <- c("low_p" = 0.15, "high_p" = 0.65)
failure_time_model <- c("correct_FG", "misspec_FG")
censoring_type <- c("none", "exponential", "curvy_uniform")

# Prediction settings
pred_timepoints <- c(0, 0.25, 0.5, 0.75, seq(1, 5, by = 0.5))
extra_settings <- list(
  tar_target(reference_patients, expand.grid(X = c(0, 1), Z = c(0, 1))),
  # We set some of the varying parameters as targets, so we can use dynamic branching later
  tar_target(failure_time_model_dyn, failure_time_model),
  tar_target(censoring_type_dyn, censoring_type)
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
        "betas" = c(0.75, 0.5),
        "base_rate" = 1,
        "base_shape" = 0.75
      )
    )
  ),
  # Generate large dataset (which we keep as a target so that we can show
  # non-param cumincs)
  tar_target(
    largedat_correct_FG,
    generate_dataset(
      n = 1e6,
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

  # FG least-false parameters depend on censoring - recover them
  tar_target(
    weibull_FG_lfps,
    recover_FG_lps(
      censoring_type = censoring_type_dyn,
      params = params_weibull_lfps,
      large_dat = generate_dataset(
        n = 1e6,
        args_event_times = list(
          mechanism = "misspec_FG",
          params = params_weibull_lfps,
          censoring_type = censoring_type_dyn
        ),
        args_missingness = list(mech_params = list("prob_missing" = 0))
      )
    ),
    pattern = map(censoring_type_dyn)
  ),

  # This are the actual simulation replications, iterate over the remaining scenario parameters
  tar_map_rep(
    name = simreps,
    values = expand.grid(
      "failure_time_model" = failure_time_model,
      "censoring_type" = censoring_type,
      stringsAsFactors = FALSE
    ),
    command = one_replication_cens_known( # change!!
      args_event_times = list(
        mechanism = failure_time_model,
        censoring_type = censoring_type,
        params = switch(
          failure_time_model,
          "correct_FG" = true_params_correct_FG,
          "misspec_FG" = params_weibull_lfps
        )
      ),
      args_missingness = list(mech_params = list("prob_missing" = 0.4, "mechanism_expr" = "Z")),
      args_imputations = list(m = 10, iters = 20, rjlimit = 1000), # m= 20
      args_predictions = list(timepoints = pred_timepoints),
      true_betas = switch(
        failure_time_model,
        "correct_FG" = true_params_correct_FG[["cause1"]][["betas"]],
        "misspec_FG" = weibull_FG_lfps[weibull_FG_lfps[["censoring_type"]] == censoring_type, ][["coefs"]]
      )
    ) |>
      cbind(prob_space = p),
    reps = 200, # number of reps is reps * batches
    batches = 2,
    combine = TRUE
  ),

  # Calculate true cumulative incidences for all reference patients
  tar_target(
    true_cuminc,
    compute_true_cuminc(
      t = pred_timepoints,
      newdat = cbind.data.frame(X = reference_patients$X, Z = reference_patients$Z),
      params = switch(
        failure_time_model_dyn,
        "correct_FG" = true_params_correct_FG,
        "misspec_FG" = params_weibull_lfps
      ),
      model_type = failure_time_model_dyn
    ) |>
      cbind(
        "failure_time_model" = failure_time_model_dyn,
        "X" = reference_patients$X,
        "Z" = reference_patients$Z,
        "prob_space" = p
      ),
    # Check against: tidyr::crossing(failure_time_model, tar_read(reference_patients))
    pattern = cross(reference_patients, failure_time_model_dyn)
  )
)


# Try also 2 extra scenarios with big betas
# Maybe misspecified + well specified with big betas? And exp censoring?
# If misspecified - need to also generate big dat..
# For now: just well specified FG?
extra_sims <- tar_rep(
  extras,
  command = one_replication_cens_known(
    args_event_times = list(
      mechanism = "correct_FG",
      censoring_type = "none",
      params = list(
        "cause1" = list(
          "formula" = ~ X + Z,
          "betas" = c(1, 1),
          "p" = 0.15,
          "base_rate" = 1,
          "base_shape" = 0.75
        ),
        "cause2" = list(
          "formula" = ~ X + Z,
          "betas" = c(1, 1),
          "base_rate" = 1,
          "base_shape" = 0.75
        )
      )
    ),
    args_missingness = list(mech_params = list("prob_missing" = 0.4, "mechanism_expr" = "Z")),
    args_imputations = list(m = 10, iters = 20, rjlimit = 1000), # m = 20
    args_predictions = list(timepoints = pred_timepoints),
    true_betas = c(1, 1)
  ) |>
    cbind(prob_space = 0.15),
  reps = 200,
  batches = 2
)


# Censoring tests
censoring_sims <- tar_map_rep(
  name = cens_sims,
  combine = TRUE,
  values = data.frame("cens_rate" = c(0.05, 0.2, 1.5, 6)),
  # rate 1.5 ~ around 50% cens, 6 leads to 75% (still like 100 events left)
  command = one_replication_cens_known(
    args_event_times = list(
      mechanism = "correct_FG",
      censoring_type = "exponential",
      params = list(
        "cause1" = list(
          "formula" = ~ X + Z,
          "betas" = c(0.75, 0.5),
          "p" = 0.15,
          "base_rate" = 1,
          "base_shape" = 0.75
        ),
        "cause2" = list(
          "formula" = ~ X + Z,
          "betas" = c(0.75, 0.5),
          "base_rate" = 1,
          "base_shape" = 0.75
        )
      ),
      censoring_params = list(
        "exponential" = cens_rate, #0.2/0.05
        "curvy_uniform" = c(0.5, 5),
        "curvyness" = 0.3
      )
    ),
    args_missingness = list(mech_params = list("prob_missing" = 0.4, "mechanism_expr" = "Z")),
    args_imputations = list(m = 10, iters = 20, rjlimit = 1000), # m = 20
    args_predictions = list(timepoints = pred_timepoints),
    true_betas = c(0.75, 0.5)
  ) |>
    cbind(prob_space = 0.15),
  reps = 50, # Let's get this quickkkkk
  batches = 8
)


# Proving a point with nonsense imps
nonsense_sims <- tar_map_rep(
  name = simplesurv_sims,
  combine = TRUE,
  values = data.frame("cens_rate" = c(0.05, 0.2, 0.5, 1.5)),
  # rate 1.5 ~ around 50% cens, 6 leads to 75% (still like 100 events left)
  command = one_replication_nonsense(
    args_event_times = list(
      mechanism = "correct_FG",
      censoring_type = "exponential",
      params = list(
        "cause1" = list(
          "formula" = ~ X + Z,
          "betas" = c(0.75, 0.5),
          "p" = 0.15,
          "base_rate" = 1,
          "base_shape" = 0.75
        ),
        "cause2" = list(
          "formula" = ~ X + Z,
          "betas" = c(0.75, 0.5),
          "base_rate" = 1,
          "base_shape" = 0.75
        )
      ),
      censoring_params = list(
        "exponential" = cens_rate, #0.2/0.05
        "curvy_uniform" = c(0.5, 5),
        "curvyness" = 0.3
      )
    ),
    args_missingness = list(mech_params = list("prob_missing" = 0.4, "mechanism_expr" = "Z")),
    args_imputations = list(m = 10, iters = 20, rjlimit = 1000), # m = 20
    args_predictions = list(timepoints = pred_timepoints),
    true_betas = c(0.75, 0.5)
  ) |>
    cbind(prob_space = 0.15),
  reps = 50, #50, # Let's get this quickkkkk
  batches = 8
)


# Here we bring together all the simulation scenarios
list(
  extra_settings,
  simulation_pipeline,
  tar_combine(
    all_simulations,
    simulation_pipeline[["simreps"]],
    command = dplyr::bind_rows(!!!.x)
  ),
  tar_combine(
    true_cuminc_all,
    simulation_pipeline[["true_cuminc"]],
    command = dplyr::bind_rows(!!!.x)
  ),
  extra_sims,
  censoring_sims,
  nonsense_sims
  # Here we pool coefficients and predictions etc.
)


# tar_read(extras) # run this one??


# Check tar meta and object sizes
# Might need to add least-false true onto the dataframe
# Calculate empirical SEs of absolute risk predictions
