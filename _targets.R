# Workhorse packages
library("targets")
library("tarchetypes")
library("future")
library("future.callr")
library("here")

# Source other packages/helper functions
source(here("packages.R"))
invisible(lapply(list.files(here("R"), full.names = TRUE), source))

# (!) Make this into targets markdown?

# To run pipeline in parallel
plan(callr)

# Data-generating parameters depend on p-space domination, so we need to iterate
# over that parameter separately
prob_space_domin <- c("low_p" = 0.15, "high_p" = 0.65)
failure_time_model <- c("correct_FG", "misspec_FG")
censoring_type <- c("none", "exponential", "curvy_uniform")

# We set some of the varying parameters as targets, so we can use dynamic branching later
dynamic_settings <- list(
  tar_target(reference_patients, expand.grid(X = c(0, 1), Z = c(0, 1))),
  tar_target(failure_time_model_dyn, failure_time_model),
  tar_target(censoring_type_dyn, censoring_type)
)

# Other global settings:
pred_timepoints <- c(0.001, 0.25, 0.5, 0.75, seq(1, 5, by = 0.5))
num_imputations <- 2#30
num_cycles <- 5 #20
num_replications <- 2 #500
num_batches <- 1
reps_per_batch <- ceiling(num_replications / num_batches)
size_data_lfps <- 1e6 # size of dataset to estimate least-false parameters
prop_missing <- 0.4


# Here we start the pipeline
simulation_pipeline_main <- tar_map(
  unlist = FALSE,
  # We iterate over the values of p
  values = list("p" = prob_space_domin),
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
  # non-param cumincs later)
  tar_target(
    largedat_correct_FG,
    generate_dataset(
      n = size_data_lfps,
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
        n = size_data_lfps,
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

  # These are the actual core simulation replications, iterate over the remaining scenario parameters
  tar_map_rep(
    name = simreps,
    values = expand.grid(
      "failure_time_model" = failure_time_model,
      "censoring_type" = censoring_type,
      stringsAsFactors = FALSE
    ),
    command = one_replication(
      args_event_times = list(
        mechanism = failure_time_model,
        censoring_type = censoring_type,
        params = switch(
          failure_time_model,
          "correct_FG" = true_params_correct_FG,
          "misspec_FG" = params_weibull_lfps
        )
      ),
      args_missingness = list(mech_params = list("prob_missing" = prop_missing, "mechanism_expr" = "1.5 * Z")),
      args_imputations = list(m = num_imputations, iters = num_cycles, rjlimit = 1000),
      args_predictions = list(timepoints = pred_timepoints),
      true_betas = switch(
        failure_time_model,
        "correct_FG" = true_params_correct_FG[["cause1"]][["betas"]],
        "misspec_FG" = weibull_FG_lfps[weibull_FG_lfps[["censoring_type"]] == censoring_type, ][["coefs"]]
      )
    ) |>
      cbind(prob_space = p),
    reps = reps_per_batch, # Total number of replications = reps * batches
    batches = num_batches,
    combine = TRUE
  ),

  # Calculate true cumulative incidences for all reference patients
  tar_target(
    true_cuminc,
    compute_true(
      t = pred_timepoints,
      newdat = reference_patients,
      params = switch(
        failure_time_model_dyn,
        "correct_FG" = true_params_correct_FG,
        "misspec_FG" = params_weibull_lfps
      ),
      model_type = failure_time_model_dyn
    ) |>
      cbind("failure_time_model" = failure_time_model_dyn, "prob_space" = p),
    pattern = map(failure_time_model_dyn)
  )
)


# We vary the censoring rate further for a subset of scenarios
# (since performance of misspecified imp model could depend on censoring)
# leads to approx 15% an 50% censored (so we have 0, 0.15, 0.3, and 0.5 as prop censored overall)
extra_cens_rates <- c(0.17, 1.44)

# simulations_censoring <- tar_map(
#   values = list("cens_rate" = extra_cens_rates),
#   unlist = FALSE,
#
#   # Same story with the least-false parameters depending on censoring,
#   # although here it will not matter so much since p = 0.65
#   tar_target(
#     weibull_FG_lfps_extracens,
#     recover_FG_lps(
#       censoring_type = "exponential",
#       params = params_weibull_lfps_0.65,
#       large_dat = generate_dataset(
#         n = size_data_lfps,
#         args_event_times = list(
#           mechanism = "misspec_FG",
#           params = params_weibull_lfps_0.65,
#           censoring_type = "exponential",
#           censoring_params = list("exponential" = cens_rate)
#         ),
#         args_missingness = list(mech_params = list("prob_missing" = 0))
#       )
#     )$coefs
#   ),
#
#   # Now run the extra replications for two scenarios, since that is where
#   # difference were most pronounce in the original replications..
#   # - p = 0.15, correct_FG
#   # - p = 0.65, misspec_FH
#   tar_map_rep(
#     name = simreps_cens,
#     combine = TRUE,
#     values = data.frame("failure_time_model" = failure_time_model),
#     command = one_replication(
#       args_event_times = list(
#         mechanism = failure_time_model,
#         censoring_type = "exponential",
#         params = switch(
#           failure_time_model,
#           "correct_FG" = true_params_correct_FG_0.15,
#           "misspec_FG" = params_weibull_lfps_0.65
#         ),
#         censoring_params = list("exponential" = cens_rate)
#       ),
#       args_missingness = list(mech_params = list("prob_missing" = prop_missing, "mechanism_expr" = "1.5 * Z")),
#       args_imputations = list(m = num_imputations, iters = num_cycles, rjlimit = 1000),
#       args_predictions = list(timepoints = pred_timepoints),
#       true_betas = switch(
#         failure_time_model,
#         "correct_FG" = true_params_correct_FG_0.15[["cause1"]][["betas"]],
#         "misspec_FG" = weibull_FG_lfps_extracens
#       )
#     ) |>
#       cbind(
#         prob_space = switch(failure_time_model, "correct_FG" = 0.15, "misspec_FG" = 0.65),
#         cens_rate = cens_rate,
#         censoring_type = "exponential"
#       ),
#     reps = reps_per_batch,
#     batches = num_batches
#   )
# )


# Additional scenarios with large covariate effects, large p; to showcase
# largest differences between smcfcs and mice. To check: use pred mean matching here?

# MCSE will be much smaller, so we can do probably do just half the reps
# stress_test <- tar_map_rep(
#   name = big_betas,
#   combine = TRUE,
#   values = data.frame("censoring_type" = c("none", "exponential")), # for both no, and exponential censoring
#   command = one_replication(
#     args_event_times = list(
#       mechanism = "correct_FG",
#       censoring_type = censoring_type,
#       censoring_params = list("exponential" = 0.25), # to get 30% cens
#       params = list(
#         "cause1" = list(
#           "formula" = ~ X + Z,
#           "betas" = c(1, 1),
#           "p" = 0.65,
#           "base_rate" = 1,
#           "base_shape" = 0.75
#         ),
#         "cause2" = list(
#           "formula" = ~ X + Z,
#           "betas" = c(1, 1),
#           "base_rate" = 1,
#           "base_shape" = 0.75
#         )
#       )
#     ),
#     args_missingness = list(mech_params = list("prob_missing" = prop_missing, "mechanism_expr" = "1.5 * Z")),
#     args_imputations = list(m = num_imputations, iters = num_cycles, rjlimit = 10000), # Large rjlimit to avoid warnings
#     args_predictions = list(timepoints = pred_timepoints),
#     args_covariates = list("X_type" = "normal"),
#     true_betas = c(1, 1)
#   ) |>
#     cbind(
#       prob_space = 0.65,
#       failure_time_model = "correct_FG"
#     ),
#   reps = 25, # 100 test replications
#   batches = 4
# )


# Here we bring together all the simulation scenarios
list(
  dynamic_settings,
  simulation_pipeline_main,
  #simulations_censoring,
  tar_combine(
    simulations_main,
    simulation_pipeline_main[["simreps"]],
    command = dplyr::bind_rows(!!!.x)
  )#,
  #tar_combine(
  #  simulations_cens,
  #  simulations_censoring[["simreps_cens"]],
  #  command = dplyr::bind_rows(!!!.x)
  #)#,
  #stress_test,
  # In reporting: forget about admin cens; just mention in-text

  # Pool predictions just for main simulations?
  # tar_target(
  #   pooled_preds_main,
  #   pool_nested_predictions(
  #     nested_df = simulations_main,
  #     new_pat = c(X = reference_patients$X, Z = reference_patients$Z)
  #   ) |>
  #     cbind(
  #       "X" = reference_patients$X,
  #       "Z" = reference_patients$Z
  #     ),
  #   pattern = map(reference_patients)
  # )
)

# To-do:
# - Calculate empirical SEs of absolute risk predictions
