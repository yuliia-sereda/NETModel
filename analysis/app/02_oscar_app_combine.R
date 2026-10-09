#!/usr/bin/env Rscript

# ==============================================================================
# App 02 - Combine the array outputs and summarise outcomes
#
# Checks that every parameter set, design point and treatment sequence is
# present, then writes:
#   app_parameter_draws.rds      the sampled parameter values of each set
#   app_sets_part_XX.rds         set-level outcomes (split for GitHub)
#   app_outcomes_summary.rds/csv mean and 95% interval across parameter sets
#                                for each design point, sequence and outcome
# Submitted by run_oscar_app_pipeline.sh after all array tasks succeed.
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(NECModel)
})

n_sets <- as.integer(Sys.getenv("NEC_N_SETS", unset = NA))
output_dir <- Sys.getenv("NEC_OUTPUT_DIR")
design_file <- Sys.getenv("NEC_DESIGN_FILE")
if (output_dir == "" || design_file == "" || is.na(n_sets)) {
  stop("Environment variables NEC_OUTPUT_DIR, NEC_DESIGN_FILE and NEC_N_SETS must be defined.")
}
output_dir <- normalizePath(output_dir, mustWork = TRUE)
design <- as.data.table(readRDS(normalizePath(design_file, mustWork = TRUE)))

part_files <- list.files(file.path(output_dir, "parts"), pattern = "^app_sets_task_\\d+\\.rds$",
                         full.names = TRUE)
message("Combining ", length(part_files), " task files")
sets <- rbindlist(lapply(part_files, readRDS), use.names = TRUE)

# Completeness: every set x design point x sequence exactly once
n_tx <- length(treatment_sequences())
expected_rows <- n_sets * nrow(design) * n_tx
if (!setequal(unique(sets$set_id), seq_len(n_sets)) || nrow(sets) != expected_rows ||
    anyDuplicated(sets, by = c("set_id", "design_id", "tx"))) {
  stop("Incomplete or duplicated results: ", nrow(sets), " rows for ",
       uniqueN(sets$set_id), " sets (expected ", expected_rows, " rows for ", n_sets, " sets).")
}

# Sampled parameter values, one row per set
param_cols <- grep("^lambda_", names(sets), value = TRUE)
draws <- unique(sets[, c("set_id", param_cols), with = FALSE])
setorder(draws, set_id)
saveRDS(draws, file.path(output_dir, "app_parameter_draws.rds"))

# Set-level outcomes (design values are in the design file), split for GitHub
outcome_cols <- app_outcome_info()$outcome
sets_out <- sets[, c("set_id", "design_id", "tx", outcome_cols), with = FALSE]
setorder(sets_out, design_id, set_id, tx)
old_parts <- list.files(output_dir, pattern = "^app_sets_part_\\d+\\.rds$", full.names = TRUE)
file.remove(old_parts)
# ~60 MB in memory per part keeps each compressed file well under GitHub's
# 100 MB limit (the outcome values are mostly unique, so they compress little)
n_parts <- ceiling(as.numeric(object.size(sets_out)) / (1024^2) / 60)
rows_per_part <- ceiling(nrow(sets_out) / n_parts)
for (i in seq_len(n_parts)) {
  idx <- seq((i - 1) * rows_per_part + 1, min(i * rows_per_part, nrow(sets_out)))
  saveRDS(sets_out[idx], file.path(output_dir, sprintf("app_sets_part_%02d.rds", i)))
}
message("Saved set-level outcomes in ", n_parts, " file(s)")

# Mean and 95% interval across parameter sets
summary <- summarise_app_outcomes(sets, design = design)
setorder(summary, design_id, tx, outcome)
saveRDS(summary, file.path(output_dir, "app_outcomes_summary.rds"))
fwrite(summary, file.path(output_dir, "app_outcomes_summary.csv"))
message("Saved outcome summary: ", nrow(summary), " rows")
message("Completed at ", Sys.time())
