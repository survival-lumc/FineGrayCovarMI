# Workhorse packages
library("targets")
library("tarchetypes")
library("future")
library("future.callr")
library("here")

# Source other packages/helper functions
source(here("packages.R"))
invisible(lapply(list.files(here("R"), full.names = TRUE), source))

tar_option_set(error = "continue")

# (!) Make this into targets markdown?

# To run pipeline in parallel
plan(callr)

# Data-generating parameters depend on p-space domination, so we need to iterate
# over that parameter separately
prob_space_domin <- c("low_p" = 0.15, "high_p" = 0.65)
failure_time_model <- c("correct_FG", "misspec_FG")
censoring_type <- c("none", "exponential", "admin")

# We set some of the varying parameters as targets, so we can use dynamic branching later
dynamic_settings <- list(
  tar_target(reference_patients, data.frame(X = c(0, 1), Z = c(0, 1))),
  tar_target(failure_time_model_dyn, failure_time_model),
  tar_target(censoring_type_dyn, censoring_type)
)

# Other global settings:
pred_timepoints <- c(0.001, 0.25, 0.5, 0.75, seq(1, 5, by = 0.5))
num_imputations <- 30
num_cycles <- 20
num_replications <- 500
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
  # Generate large dataset + estimate least-false Weibull parameters, for use in the misspecified FG scenario
  tar_target(
    params_weibull_lfps,
    recover_weibull_lfps(
      large_dat = generate_dataset(
        n = size_data_lfps,
        args_event_times = list(
          mechanism = "correct_FG",
          params = true_params_correct_FG,
          censoring_type = "none"
        ),
        args_missingness = list(mech_params = list("prob_missing" = 0))
      ),
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
      args_imputations = list(
        m = num_imputations,
        iters = num_cycles,
        rjlimit = 1000,
        rhs_kmi = "1"
      ),
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

# Summarize the simulations
summarized_sims <- list(
  dynamic_settings,
  simulation_pipeline_main,
  tar_combine(true_cuminc_all, simulation_pipeline_main[["true_cuminc"]]),
  tar_combine(
    simulations_main,
    simulation_pipeline_main[["simreps"]],
    command = {
      setDT(vctrs::vec_c(!!!.x))[, .(
        coefs = list(rbindlist(coefs_summary, idcol = "rep_id")),
        preds = list(rbindlist(preds_summary, idcol = "rep_id"))
      ), by = c("prob_space","failure_time_model", "censoring_type", "method")]
    }
  ),
  tar_target(
    coefs_main,
    command = {
      df <- rbindlist(
        with(
          simulations_main,
          Map(
            f = cbind,
            method = method,
            coefs,
            prob_space = prob_space,
            failure_time_model = failure_time_model,
            censoring_type = censoring_type
          )
        ), fill = TRUE
      )
      df[, term := ifelse(grepl(pattern = "^X", term), "X", as.character(term))]
    }
  ),
  tar_target(
    preds_main,
    command = {
      df <- rbindlist(
        with(
          simulations_main,
          Map(
            f = cbind,
            method = method,
            preds,
            prob_space = prob_space,
            failure_time_model = failure_time_model,
            censoring_type = censoring_type
          )
        ), fill = TRUE
      )
      df[is.na(imp), imp := 0]
    }
  ),
  tar_target(
    pooled_preds_main,
    pool_nested_predictions(preds_main, new_pat = as.numeric(reference_patients)),
    pattern = map(reference_patients)
  )
)

applied_imp_settings <- list(
  num_imputations = 100,
  num_cycles = 20,
  num_batches = 10,
  rjlimit = 10000,
  rhs_cens = "year_allo1_decades",
  cause = 1 # relapse
)

# Reminder total number is batches/reps
applied_example <- list(
  tar_target(file, "data-raw//dat_clean.rds", format = "file"),
  tar_target(applied_dat_raw, data.table(readRDS(file))),
  tar_target(applied_dat, process_applied_dat(applied_dat_raw)),
  tar_rep(
    applied_impdats,
    one_imputation_applied_dat(dat_processed = applied_dat, imp_settings = applied_imp_settings),
    reps = ceiling(applied_imp_settings$num_imputations / applied_imp_settings$num_batches),
    batches = applied_imp_settings$num_batches, # for parallelizing
    format = "fst"
  ),
  tar_target(applied_dat_pooled, pool_applied_dat(applied_impdats, applied_dat))
)


# Here we bring together all the simulation scenarios
list(
  applied_example,
  summarized_sims,
  tar_quarto(simulation_results, path = "analysis/simulation-results.qmd"),
  tar_quarto(supplement, path = "analysis/supplementary-material.qmd")
)
