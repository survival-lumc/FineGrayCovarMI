# Workhorse packages
library("targets")
library("tarchetypes")
library("future")
library("future.callr")
library("here")

# Source the rest
source(here("packages.R"))
source(here("R/data-generation-new.R")) # or
#functions <- list.files(here("R"), full.names = TRUE)
#invisible(lapply(functions, source))

# MAKE INTO TARGETS MARKDOWN!

# See https://books.ropensci.org/targets/hpc.html#future for parallelizing
plan(callr)

#tar_option_set(packages = project_pkgs)



# Start of pipeline
list(

)
# list(
#   tar_target(data, data.frame(x = sample.int(100), y = sample.int(100))),
#   tar_target(summary, summ(data)) # Call your custom functions as needed.
# )
