# App 04 - Prepare the data for the clinician app
#
# Run from the repository root after the app simulations
# (run_oscar_app_pipeline.sh) have been copied to Data/app_outputs/COMPETE/.
# Writes small, app-ready files to app/data/, so the app depends only on CRAN
# packages and not on NECModel or its GitHub dependencies.

devtools::load_all()
library(data.table)

input_dir <- file.path("Data", "app_outputs", "COMPETE")
out_dir <- file.path("app", "data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Design grid, in the units shown in the app
# ------------------------------------------------------------------------------
design <- as.data.table(readRDS(file.path("design", "app_design_grid.rds")))
design[, `:=`(
  loss2 = as.integer(round(step2_attenuation * 100)),
  loss3 = as.integer(round(step3_attenuation * 100)),
  mds = aml_mds_multiplier
)]
saveRDS(design[, .(design_id, loss2, loss3, mds)], file.path(out_dir, "design.rds"))

# ------------------------------------------------------------------------------
# Outcome definitions
# ------------------------------------------------------------------------------
info <- app_outcome_info()
info[, group := fcase(startsWith(outcome, "os_"), "Overall survival",
                      startsWith(outcome, "pfs_"), "Progression-free survival",
                      default = "MDS/AML")]
saveRDS(info, file.path(out_dir, "outcome_info.rds"))

# ------------------------------------------------------------------------------
# Mean and 95% range across parameter sets
# ------------------------------------------------------------------------------
outcomes <- readRDS(file.path(input_dir, "app_outcomes_summary.rds"))
saveRDS(outcomes[, .(design_id, tx, outcome, mean, lower95, upper95)], file.path(out_dir, "outcomes.rds"))

# ------------------------------------------------------------------------------
# Differences from every possible reference strategy
# ------------------------------------------------------------------------------
differences <- readRDS(file.path(input_dir, "app_differences_summary.rds"))
forward <- differences[, .(design_id, reference = tx2, tx = tx1, outcome, diff_unit,
                           mean_diff, lower95, upper95, p_better = p_tx1_better)]
# The reverse direction: tx2 - tx1 = -(tx1 - tx2)
reverse <- differences[, .(design_id, reference = tx1, tx = tx2, outcome, diff_unit,
                           mean_diff = -mean_diff, lower95 = -upper95, upper95 = -lower95,
                           p_better = 1 - p_tx1_better)]
saveRDS(rbind(forward, reverse), file.path(out_dir, "vs_reference.rds"))

# ------------------------------------------------------------------------------
# Set-level results for the uncertainty regions
# ------------------------------------------------------------------------------
part_files <- list.files(input_dir, pattern = "^app_sets_part_\\d+\\.rds$", full.names = TRUE)
sets <- rbindlist(lapply(part_files, readRDS))
stopifnot(nrow(sets) == nrow(design) * 1000 * length(treatment_sequences()))

# Mean OS, mean PFS and 5-year MDS/AML CIR per set, for the uncertainty regions
clouds <- sets[, .(set_id = as.integer(set_id), design_id = as.integer(design_id), tx = factor(tx),
                   os_mean = round(os_mean, 4), pfs_mean = round(pfs_mean, 4),
                   mds_aml_5y_cir = round(mds_aml_5y_cir, 5))]
setkey(clouds, design_id)
saveRDS(clouds, file.path(out_dir, "clouds.rds"), compress = "xz")

# ------------------------------------------------------------------------------
# Model inputs for the methods page
# ------------------------------------------------------------------------------
param_input <- fread(file.path("Data", "input_params_COMPETE.csv"))
param_labels <- c(
  lambda_everolimus_12 = "Progression on first-line everolimus",
  lambda_prrt_12 = "Progression on first-line PRRT",
  lambda_cabozantinib_23 = "Progression on second-line cabozantinib",
  lambda_46 = "Death after progression on salvage treatment",
  lambda_dot5 = "MDS/AML per PRRT dose",
  lambda_57 = "Death after MDS/AML"
)
inputs <- param_input[param %in% names(param_labels), .(
  transition = param_labels[param],
  estimate = fifelse(point_stat == "median",
                     sprintf("Median %.1f months", est * 12),
                     sprintf("%.4f per year per dose", est)),
  source
)]
inputs[, transition := factor(transition, levels = param_labels)]
setorder(inputs, transition)
inputs[, transition := as.character(transition)]
saveRDS(inputs, file.path(out_dir, "inputs.rds"))

message("App data written to ", out_dir, ":")
print(data.table(file = list.files(out_dir), size_MB = round(file.size(list.files(out_dir, full.names = TRUE)) / 1024^2, 2)))
