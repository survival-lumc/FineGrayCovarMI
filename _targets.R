# Workhorse packages
library(targets)
library(tarchetypes)
library(future)
library(future.callr)

# See https://books.ropensci.org/targets/hpc.html#future for parallelizing
plan(callr)

# All packages used by the projects
project_pkgs <- c(
  "data.table",
  "prodlim",
  "riskRegression",
  "survival",
  "mice",
  "ggplot2",
  "kableExtra"
)

tar_option_set(packages = project_pkgs)
# Uncomment if running scripts interactively:
# sapply(project_pkgs, function(pkg) require(pkg, character.only = TRUE))


# Start pipeline here:
# list(
#   tar_target(data, data.frame(x = sample.int(100), y = sample.int(100))),
#   tar_target(summary, summ(data)) # Call your custom functions as needed.
# )
