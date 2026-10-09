#!/usr/bin/env Rscript

# ==============================================================================
# App 01 - Simulate the app design grid under parameter uncertainty (array task)
#
# Draws n_sets sets of model parameters (the same draws in every task) and, for
# this task's share of the sets, simulates the 100,000-person population at
# every point of the app design grid. Outcomes per set: mean and median OS and
# PFS, OS and PFS at 1/3/5 years, MDS/AML cumulative incidence at 1/3/5 years.
# Submitted as a SLURM array by run_oscar_app_pipeline.sh.
# Output: <NEC_OUTPUT_DIR>/parts/app_sets_task_XX.rds
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(NECModel)
})

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA))
n_tasks <- as.integer(Sys.getenv("NEC_N_TASKS", unset = NA))
n_sets <- as.integer(Sys.getenv("NEC_N_SETS", unset = NA))
seed <- as.integer(Sys.getenv("NEC_SEED", unset = NA))

paths <- c(
  NEC_INPUT_FILE  = Sys.getenv("NEC_INPUT_FILE"),
  NEC_DESIGN_FILE = Sys.getenv("NEC_DESIGN_FILE"),
  NEC_OUTPUT_DIR  = Sys.getenv("NEC_OUTPUT_DIR")
)
if (any(paths == "") || anyNA(c(task_id, n_tasks, n_sets, seed))) {
  stop("Environment variables NEC_INPUT_FILE, NEC_DESIGN_FILE, NEC_OUTPUT_DIR, ",
       "NEC_N_TASKS, NEC_N_SETS, NEC_SEED and SLURM_ARRAY_TASK_ID must be defined.")
}
input_file  <- normalizePath(paths[["NEC_INPUT_FILE"]], mustWork = TRUE)
design_file <- normalizePath(paths[["NEC_DESIGN_FILE"]], mustWork = TRUE)
parts_dir   <- file.path(normalizePath(paths[["NEC_OUTPUT_DIR"]], mustWork = TRUE), "parts")
dir.create(parts_dir, showWarnings = FALSE)

# This task's share of the parameter sets
if (n_tasks > n_sets || task_id < 1 || task_id > n_tasks) {
  stop("Need 1 <= task_id <= n_tasks <= n_sets (task ", task_id, ", ", n_tasks, " tasks, ", n_sets, " sets).")
}
set_ids <- split(seq_len(n_sets), cut(seq_len(n_sets), n_tasks, labels = FALSE))[[task_id]]

message("Task ", task_id, " of ", n_tasks, ": parameter sets ", min(set_ids), "-", max(set_ids),
        " of ", n_sets, " on ", n_cores, " cores")
message("Input parameters: ", input_file)
message("Design file:      ", design_file)

pop <- make_population(seed = 5678)
param_input <- as.data.table(read.csv(input_file))
design <- as.data.table(readRDS(design_file))
message("Design points: ", nrow(design))
message("Started at ", Sys.time())

result <- sim_parameter_sets(
  param_input = param_input,
  design = design[, .(discount1, discount2, aml_mds_multiplier)],
  population = pop,
  treatments = treatment_sequences(),
  n_sets = n_sets,
  seed = seed,
  n_cores = n_cores,
  set_ids = set_ids,
  summarise = summarise_population_app
)

output_file <- file.path(parts_dir, sprintf("app_sets_task_%02d.rds", task_id))
saveRDS(result, file = output_file)
message("Saved ", nrow(result), " rows to ", output_file, " at ", Sys.time())
