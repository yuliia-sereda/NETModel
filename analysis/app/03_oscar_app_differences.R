#!/usr/bin/env Rscript

# ==============================================================================
# App 03 - Pairwise differences between treatment sequences (dominance analysis)
#
# For each design point, pair of sequences (tx1, tx2) and outcome: the
# difference tx1 - tx2 within each parameter set, summarised by its mean, 95%
# interval and the share of sets in which tx1 is better.
# Submitted by run_oscar_app_pipeline.sh after 02 succeeds.
# Output: <NEC_OUTPUT_DIR>/app_differences_summary.rds/csv
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(NECModel)
})

data.table::setDTthreads(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1")))

output_dir <- Sys.getenv("NEC_OUTPUT_DIR")
design_file <- Sys.getenv("NEC_DESIGN_FILE")
if (output_dir == "" || design_file == "") {
  stop("Environment variables NEC_OUTPUT_DIR and NEC_DESIGN_FILE must be defined.")
}
output_dir <- normalizePath(output_dir, mustWork = TRUE)
design <- as.data.table(readRDS(normalizePath(design_file, mustWork = TRUE)))

part_files <- list.files(output_dir, pattern = "^app_sets_part_\\d+\\.rds$", full.names = TRUE)
if (length(part_files) == 0L) {
  stop("No set-level outcome files (app_sets_part_XX.rds) in ", output_dir)
}
sets <- rbindlist(lapply(part_files, readRDS), use.names = TRUE)
message("Set-level rows: ", nrow(sets))

differences <- summarise_app_differences(sets, design = design)
setorder(differences, design_id, tx1, tx2, outcome)
saveRDS(differences, file.path(output_dir, "app_differences_summary.rds"))
fwrite(differences, file.path(output_dir, "app_differences_summary.csv"))
message("Saved differences summary: ", nrow(differences), " rows")
message("Completed at ", Sys.time())
