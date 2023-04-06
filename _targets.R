# Workhorse packages
library("targets")
library("tarchetypes")
library("future")
library("future.callr")
library("here")

# Source the rest
source(here("packages.R"))
source(here("R/data-generation.R"))
source(here("R/smcfcs.finegray.R"))
source(here("R/crprep-timefixed.R"))
# or
#functions <- list.files(here("R"), full.names = TRUE)
#invisible(lapply(functions, source))

# MAKE INTO TARGETS MARKDOWN!

# See https://books.ropensci.org/targets/hpc.html#future for parallelizing
plan(callr)

#tar_option_set(packages = project_pkgs)

# use variance = FALSE in FGR to make it faster?


# Or maybe only part for the LFPs (without censoring bit)
scenarios <- expand.grid(
  "failure_time_model" = c("correct_FG", "misspec_FG"),
  "prob_space_domin" = c("low_p" = 0.25, "high_p" = 0.75),
  "censoring_type" = c("none", "exponential", "uniform")
)



# Start of pipeline, iterate over scenarios df
# list(
#   tar_map(
#     value = list("prob_space_domin" = c("low_p" = 0.25, "high_p" =0.75)),
#     tar_target(
#       true_params_correct_FG,
#       list(
#         "cause1" = list(
#           "formula" = ~ X + Z,
#           "betas" = c(0.75, 0.5),
#           "p" = prob_space_domin,
#           "base_rate" = 1,
#           "base_shape" = 0.75
#         ),
#         "cause2" = list(
#           "formula" = ~ X + Z,
#           "betas" = c(0.75, 0.5),
#           "base_rate" = 1,
#           "base_shape" = 0.75
#         )
#       )
#     )
#   )
# )

# First get it working for one p
list(
  tar_target(
    true_params_correct_FG,
    list(
      "cause1" = list(
        "formula" = ~ X + Z,
        "betas" = c(0.75, 0.5),
        "p" = 0.25,
        "base_rate" = 1,
        "base_shape" = 0.75
      ),
      "cause2" = list(
        "formula" = ~ X + Z,
        "betas" = c(0.75, 0.5), # pick some that result in non-zero cause-spec coefs!!
        "base_rate" = 1, # might need to increase this? nah probs wont matter
        "base_shape" = 0.75
      )
    )
  ),
  # Need to do this for large and small "p" (domination of prob space)
  # Can probably to some sort of apply to vector of p
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
  tar_target(
    params_weibull_lfps,
    recover_weibull_lfps(
      large_dat = largedat_correct_FG,
      params_correct_FG = true_params_correct_FG
    )
  ),
  tar_target(
    largedat_weibull_cause_spec,
    generate_dataset(
      n = 1e6,
      args_event_times = list(
        mechanism = "misspec_FG",
        params = params_weibull_lfps,
        censoring_type = "none"
      ),
      args_missingness = list(mech_params = list("prob_missing" = 0))
    )
  ),
  tar_target(
    largedat_weibull_FG_lfps,
    recover_fg_lps(largedat_weibull_cause_spec)
  )
  # 'True' least false in cause-spec.. should depend on censoring?
  # Also make scenarios table..
)

