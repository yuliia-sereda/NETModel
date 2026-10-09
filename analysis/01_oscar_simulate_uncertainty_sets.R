#!/usr/bin/env Rscript

# ==============================================================================
# 01 Simulate NEC data on OSCAR - Uncertainty Sets
#
# Simulates all treatment sequences at every point of the maxPro design.
# Submitted by run_oscar_pipeline.sh; paths are passed as environment variables.
# Output: <NEC_OUTPUT_DIR>/NEC_data_uncertainty_sets.rds (input to 02)
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(NECModel)
})

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
message("Cores: ", n_cores)

# Read paths passed from the SLURM script
paths <- c(
  NEC_INPUT_FILE  = Sys.getenv("NEC_INPUT_FILE"),
  NEC_DESIGN_FILE = Sys.getenv("NEC_DESIGN_FILE"),
  NEC_OUTPUT_DIR  = Sys.getenv("NEC_OUTPUT_DIR")
)
if (any(paths == "")) {
  stop("Environment variable(s) not defined: ", paste(names(paths)[paths == ""], collapse = ", "))
}
input_file  <- normalizePath(paths[["NEC_INPUT_FILE"]], mustWork = TRUE)
design_file <- normalizePath(paths[["NEC_DESIGN_FILE"]], mustWork = TRUE)
output_dir  <- normalizePath(paths[["NEC_OUTPUT_DIR"]], mustWork = TRUE)

message("Input parameters: ", input_file)
message("Design file:      ", design_file)
message("Output directory: ", output_dir)

# ------------------------------------------------------------------------------
# Model inputs
# ------------------------------------------------------------------------------

seed <- 5678
pop <- make_population(seed = seed)
param_input <- as.data.table(read.csv(input_file))
D <- as.data.table(readRDS(design_file))

message("Number of design rows: ", nrow(D))
message("Starting uncertainty-set simulation at ", Sys.time())

# ------------------------------------------------------------------------------
# Simulate uncertainty sets
# ------------------------------------------------------------------------------

NEC_data <- sim_data(
  param_input = param_input,
  design = D,
  population = pop,
  treatments = treatment_sequences(),
  seed = seed,
  show_progress = FALSE,
  n_cores = n_cores
)

output_file <- file.path(output_dir, "NEC_data_uncertainty_sets.rds")
saveRDS(NEC_data, file = output_file)

message("Saved uncertainty-set results to: ", output_file)
message("Simulation completed at ", Sys.time())

sessionInfo()
