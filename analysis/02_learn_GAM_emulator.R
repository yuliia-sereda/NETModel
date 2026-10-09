#!/usr/bin/env Rscript

# ==============================================================================
# 02 Learn GAM emulators of the uncertainty-set simulation
#
# Fits one GAM per treatment sequence and outcome over the design parameters
# and predicts on a regular grid of the design space. Fits run in parallel.
# Submitted by run_oscar_pipeline.sh after 01 (sets) succeeds.
# Input:  <NEC_OUTPUT_DIR>/NEC_data_uncertainty_sets.rds
# Output: <NEC_OUTPUT_DIR>/NEC_GAM_predictions.rds (input to 03)
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(mgcv)
})

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
message("Cores: ", n_cores)

output_dir <- Sys.getenv("NEC_OUTPUT_DIR")
if (output_dir == "") {
  stop("Environment variable NEC_OUTPUT_DIR is not defined.")
}
output_dir <- normalizePath(output_dir, mustWork = TRUE)
input_file <- normalizePath(file.path(output_dir, "NEC_data_uncertainty_sets.rds"), mustWork = TRUE)

message("Input:            ", input_file)
message("Output directory: ", output_dir)

dat <- readRDS(input_file)

# Data for predictions: a regular grid over the uncertainty set. Each subsequent
# position can lose up to 50% of effectiveness (progression rate up to 1.5 times
# that of the previous position): 1 <= discount1 <= 1.5 and
# discount1 <= discount2 <= 1.5 * discount1.
discount1 <- seq(1.00, 1.50, by = 0.01)
discount2 <- seq(1.00, 2.25, by = 0.01)
aml_mds_multiplier <- seq(1.00, 3.00, by = 0.05)
tol <- 1e-9
newdat <- CJ(discount1, discount2, aml_mds_multiplier)[
  discount2 >= discount1 - tol & discount2 <= 1.5 * discount1 + tol
]
message("Prediction grid: ", nrow(newdat), " points")

# GAM: main effects plus ti() interactions, so each interaction term
# excludes the main effects already in the model.
uncertainty_vars <- names(newdat)
outcome_vars <- c("os_mean", "pfs_mean", "mds_aml_5y_cir")
treatments <- unique(dat$tx)
fits <- CJ(treatment = treatments, outcome = outcome_vars, sorted = FALSE)

fit_gam <- function(tx_i, o) {
  fml <- as.formula(paste0(
    o,
    " ~ s(discount1) + s(discount2) + s(aml_mds_multiplier) + ",
    "ti(discount1, discount2) + ti(discount1, aml_mds_multiplier) + ti(discount2, aml_mds_multiplier) + ",
    "ti(discount1, discount2, aml_mds_multiplier)"
  ))
  model <- bam(fml, data = dat[tx == tx_i], discrete = TRUE, nthreads = 1)
  list(
    residual_sd = sqrt(model$sig2),
    dev_expl = summary(model)$dev.expl,
    prediction = as.vector(predict(model, newdata = newdat))
  )
}

message("Fitting ", nrow(fits), " GAMs at ", Sys.time())
results <- parallel::mclapply(
  seq_len(nrow(fits)),
  function(i) fit_gam(fits$treatment[i], fits$outcome[i]),
  mc.cores = n_cores,
  mc.preschedule = FALSE
)

# A worker that errors returns a try-error; one that is killed (e.g. out of
# memory) returns NULL. Treat both as failures rather than dropping rows.
failed <- vapply(results, function(x) is.null(x) || inherits(x, "try-error"), logical(1))
if (any(failed)) {
  stop("GAM fit(s) failed for: ",
       paste(fits[failed, paste(treatment, outcome)], collapse = ", "), ": ",
       results[[which(failed)[1]]])
}

# Collect fit diagnostics and predictions (one row per grid point and treatment)
residual_sd <- copy(fits)[, `:=`(
  residual_sd = vapply(results, `[[`, numeric(1), "residual_sd"),
  dev_expl = vapply(results, `[[`, numeric(1), "dev_expl")
)]
print(residual_sd)

predictions <- rbindlist(lapply(treatments, function(tx_i) {
  pred_tx <- cbind(newdat, treatment = tx_i)
  for (o in outcome_vars) {
    pred_tx[, (o) := results[[fits[, which(treatment == tx_i & outcome == o)]]]$prediction]
  }
  pred_tx
}))

gam_results <- list(residual_sd = residual_sd, predictions = predictions)
output_file <- file.path(output_dir, "NEC_GAM_predictions.rds")
saveRDS(gam_results, file = output_file)

message("Saved GAM predictions to: ", output_file)
message("Completed at ", Sys.time())

sessionInfo()
