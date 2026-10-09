#!/usr/bin/env Rscript

# ==============================================================================
# 01 Simulate NEC data on OSCAR - Uncertainty Bounds
#
# Probabilistic sensitivity analysis of the two extreme scenarios (no
# attenuation / baseline MDS/AML risk and maximum attenuation / 3x MDS/AML
# risk), for the Pareto front: draws 1000 sets of model parameters from their
# uncertainty distributions, simulates the 100,000-person population under each
# set (everyone shares the set's parameter values), and summarises mean OS,
# mean PFS and 5-year MDS/AML CIR per set. Sets are simulated in parallel.
# Submitted by run_oscar_pipeline.sh; paths are passed as environment variables.
# Output: <NEC_OUTPUT_DIR>/NEC_data_uncertainty_bounds.rds
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(NECModel)
})

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
message("Cores: ", n_cores)

# Read paths passed from the SLURM script
paths <- c(
  NEC_INPUT_FILE = Sys.getenv("NEC_INPUT_FILE"),
  NEC_OUTPUT_DIR = Sys.getenv("NEC_OUTPUT_DIR")
)
if (any(paths == "")) {
  stop("Environment variable(s) not defined: ", paste(names(paths)[paths == ""], collapse = ", "))
}
input_file <- normalizePath(paths[["NEC_INPUT_FILE"]], mustWork = TRUE)
output_dir <- normalizePath(paths[["NEC_OUTPUT_DIR"]], mustWork = TRUE)

message("Input parameters: ", input_file)
message("Output directory: ", output_dir)

# ------------------------------------------------------------------------------
# Model inputs
# ------------------------------------------------------------------------------

pop <- make_population(seed = 5678)
param_input <- as.data.table(read.csv(input_file))
D_bounds <- data.table(
  discount1 = c(1, 1.5),
  discount2 = c(1, 2.25),
  aml_mds_multiplier = c(1, 3)
)

# ------------------------------------------------------------------------------
# Simulate uncertainty-bound data for the Pareto front
# ------------------------------------------------------------------------------

n_sets <- 1000L
message("Starting ", n_sets, " parameter sets at ", Sys.time())

# Each set has its own random-number stream; within a set both scenarios share
# random numbers.
NEC_bounds_data <- sim_parameter_sets(
  param_input = param_input,
  design = D_bounds,
  population = pop,
  treatments = treatment_sequences(),
  n_sets = n_sets,
  seed = 5678,
  n_cores = n_cores
)

output_file <- file.path(output_dir, "NEC_data_uncertainty_bounds.rds")
saveRDS(NEC_bounds_data, file = output_file)

message("Saved uncertainty-bound results to: ", output_file)
message("Simulation completed at ", Sys.time())

sessionInfo()
