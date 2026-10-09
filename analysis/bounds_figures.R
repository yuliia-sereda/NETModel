# Bounds figures
#
# Range (min-max) of the GAM-emulated mean OS, mean PFS, and 5-year MDS/AML CIR
# for each treatment sequence across the whole uncertainty set.
# Input:  Data/simulation_outputs/<set>/NEC_GAM_predictions.rds
#         (from analysis/02_learn_GAM_emulator.R)
# Output: Fig/<set>/Bounds_os_mean.pdf, Bounds_pfs_mean.pdf,
#         mean_outcomes_ranges_all_scenarios.csv

devtools::load_all()
library(data.table)
library(ggplot2)

# Parameter sets to plot; sets without GAM outputs are skipped
param_sets <- c("NETTER&RADIANT", "COMPETE")
input_root <- file.path("Data", "simulation_outputs")
fig_root <- "Fig"

outcome_vars <- c("os_mean", "pfs_mean", "mds_aml_5y_cir")
axis_labels <- c(os_mean = "Mean OS, years", pfs_mean = "Mean PFS, years")

format_range <- function(lower, upper, digits = 2) {
  paste0(round(lower, digits), "-", round(upper, digits))
}

for (param_set in param_sets) {
  input_file <- file.path(input_root, param_set, "NEC_GAM_predictions.rds")
  if (!file.exists(input_file)) {
    message("Skipping ", param_set, ": ", input_file, " not found")
    next
  }
  message("=== Parameter set: ", param_set)
  fig_dir <- file.path(fig_root, param_set)
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  dat <- readRDS(input_file)[["predictions"]]
  bounds <- get_bounds(predictions = dat, out_vars = outcome_vars, by_vars = "treatment")
  bounds[, treatment := factor(treatment, levels = names(treatment_sequences()))]
  setorder(bounds, treatment)
  print(bounds[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])

  mean_outcomes_ranges_all_scenarios <- bounds[, .(
    treatment,
    os_mean_range = format_range(os_mean..min, os_mean..max),
    pfs_mean_range = format_range(pfs_mean..min, pfs_mean..max),
    mds_aml_5y_cir_pct_range = format_range(mds_aml_5y_cir..min * 100, mds_aml_5y_cir..max * 100)
  )]
  write.csv(mean_outcomes_ranges_all_scenarios,
            file.path(fig_dir, "mean_outcomes_ranges_all_scenarios.csv"),
            row.names = FALSE)

  for (o in names(axis_labels)) {
    p <- ggplot(
      data = bounds,
      aes(
        y = treatment,
        x = .data[[paste0(o, "..min")]],
        xend = .data[[paste0(o, "..max")]]
      )
    ) +
      geom_segment(aes(yend = treatment), linewidth = 1.2, color = "#440154FF") +
      labs(title = "", y = "", x = axis_labels[[o]]) +
      # Axis fitted to this outcome's bounds, with half-year breaks
      scale_x_continuous(breaks = scales::breaks_width(0.5)) +
      theme_bw()
    ggsave(file.path(fig_dir, paste0("Bounds_", o, ".pdf")), plot = p, width = 4, height = 3)
  }
}
