#!/usr/bin/env Rscript

# ==============================================================================
# 03 Calculate pairwise treatment differences for the Hasse graphs
#
# For every pair of treatment sequences and every grid point, computes the
# difference (tx1 - tx2) in each emulated outcome.
# Submitted by run_oscar_pipeline.sh after 02 succeeds.
# Input:  <NEC_OUTPUT_DIR>/NEC_GAM_predictions.rds
# Output: <NEC_OUTPUT_DIR>/NEC_GAM_differences_part_XX.rds
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

data.table::setDTthreads(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1")))

output_dir <- Sys.getenv("NEC_OUTPUT_DIR")
if (output_dir == "") {
  stop("Environment variable NEC_OUTPUT_DIR is not defined.")
}
output_dir <- normalizePath(output_dir, mustWork = TRUE)
input_file <- normalizePath(file.path(output_dir, "NEC_GAM_predictions.rds"), mustWork = TRUE)

message("Input:            ", input_file)
message("Output directory: ", output_dir)

dat <- readRDS(input_file)[["predictions"]]

keys <- c("discount1", "discount2", "aml_mds_multiplier")
outcome_vars <- c("os_mean", "pfs_mean", "mds_aml_5y_cir")

pairs_list <- combn(unique(dat$treatment), m = 2, simplify = FALSE)

# One wide table per outcome: one row per grid point, one column per treatment
setkeyv(dat, keys)
wide <- lapply(setNames(outcome_vars, outcome_vars), function(o) {
  dcast(dat, as.formula(paste(paste(keys, collapse = " + "), "~ treatment")), value.var = o)
})
grid <- wide[[1]][, ..keys]
for (o in outcome_vars) {
  stopifnot(identical(wide[[o]][, ..keys], grid))
}

delta <- rbindlist(lapply(pairs_list, function(x) {
  out <- copy(grid)[, `:=`(tx1 = x[1], tx2 = x[2])]
  for (o in outcome_vars) {
    out[, paste0(o, "_delta") := wide[[o]][[x[1]]] - wide[[o]][[x[2]]]]
  }
  out
}))

message("Rows: ", nrow(delta))

# Split into parts small enough for GitHub (< 100 MB per file)
total_size <- as.numeric(object.size(delta)) / (1024^2)
target_size <- 200
n_parts <- ceiling(total_size / target_size)
rows_per_part <- ceiling(nrow(delta) / n_parts)
old_parts <- list.files(output_dir, pattern = "^NEC_GAM_differences_part_\\d+\\.rds$", full.names = TRUE)
file.remove(old_parts)
for (i in seq_len(n_parts)) {
  idx_start <- ((i - 1) * rows_per_part) + 1
  idx_end <- min(i * rows_per_part, nrow(delta))
  saveRDS(delta[idx_start:idx_end], file = file.path(output_dir, sprintf("NEC_GAM_differences_part_%02d.rds", i)))
}

message("Saved ", n_parts, " difference file(s) to: ", output_dir)
message("Completed at ", Sys.time())
