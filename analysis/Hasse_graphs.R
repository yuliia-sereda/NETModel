# Hasse graphs
#
# Hasse diagrams of the treatment preferences under interval dominance,
# maximality, optimism and pessimism (without and with a 2-month margin), and
# under the Hurwicz criterion over a grid of alpha values.
# Input:  Data/simulation_outputs/<set>/NEC_GAM_predictions.rds (from 02)
#         Data/simulation_outputs/<set>/NEC_GAM_differences_part_XX.rds (from 03)
# Output: Fig/<set>/Hasse_graphs_no_margin.pdf, Hasse_graphs_2mo_margin.pdf,
#         Hurwicz_graphs.pdf

devtools::load_all()
library(data.table)
library(igraph)

# Parameter sets to plot; sets without GAM outputs are skipped
param_sets <- c("NETTER&RADIANT", "COMPETE")
input_root <- file.path("Data", "simulation_outputs")
fig_root <- "Fig"

outcome_vars <- c("os_mean", "pfs_mean")
outcome_labels <- c(os_mean = "mean OS", pfs_mean = "mean PFS", mds_aml_5y_cir = "5-year MDS/AML CIR")
# Outcomes where lower is better
reversed_outcomes <- "mds_aml_5y_cir"

criteria <- list(
  "Interval dominance" = pairwise_relation_interval_dominance,
  "Maximality" = NULL, # uses the pairwise differences, not the bounds
  "Optimism" = pairwise_relation_optimism,
  "Pessimism" = pairwise_relation_pessimism
)

#' Adjacency matrix of treatment preferences for one criterion and outcome
get_adjacency <- function(criterion, outcome, margin, bounds, diffs, alpha = NULL) {
  reverse <- outcome %in% reversed_outcomes
  if (criterion == "Maximality") {
    adjacency_matrix_maximality(
      differences = diffs,
      out_var = outcome,
      digits = 3,
      reverse = reverse,
      added_margin = margin
    )
  } else {
    pw_relation <- if (criterion == "Hurwicz") pairwise_relation_hurwicz else criteria[[criterion]]
    adjacency_matrix_from_bounds(
      bounds = bounds,
      pw_relation = pw_relation,
      out_var = outcome,
      digits = 3,
      added_margin = margin,
      reverse = reverse,
      alpha = alpha
    )
  }
}

plot_hasse <- function(adjmat, main, vertex_size = 27, label_cex = 1, ...) {
  g <- hasse_graph(adjmat)
  plot(
    g,
    layout = hasse_layout(g, ...),
    main = main,
    vertex.size = vertex_size,
    vertex.label.cex = label_cex,
    vertex.label.color = "grey25",
    vertex.color = "white",
    vertex.frame.color = "grey25",
    edge.color = "grey25",
    edge.width = 1,
    edge.arrow.size = 0.5
  )
}

margin_label <- function(margin) {
  if (margin == 0) "no margin" else paste0(round(margin * 12, 1), "-month margin")
}

for (param_set in param_sets) {
  input_dir <- file.path(input_root, param_set)
  predictions_file <- file.path(input_dir, "NEC_GAM_predictions.rds")
  difference_files <- list.files(input_dir, pattern = "^NEC_GAM_differences_part_\\d+\\.rds$", full.names = TRUE)
  if (!file.exists(predictions_file) || length(difference_files) == 0L) {
    message("Skipping ", param_set, ": GAM predictions or differences not found in ", input_dir)
    next
  }
  message("=== Parameter set: ", param_set, " (", length(difference_files), " difference files)")
  fig_dir <- file.path(fig_root, param_set)
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  # For interval dominance, optimism, pessimism and Hurwicz we need the bounds;
  # for maximality, the pairwise differences
  predictions <- readRDS(predictions_file)[["predictions"]]
  bounds <- get_bounds(predictions = predictions, out_vars = outcome_vars)
  diffs <- rbindlist(lapply(difference_files, readRDS), use.names = TRUE)

  for (margin in c(0, 2 / 12)) {
    file_suffix <- if (margin == 0) "no_margin" else "2mo_margin"
    pdf(file.path(fig_dir, paste0("Hasse_graphs_", file_suffix, ".pdf")), width = 7, height = 5)
    for (criterion in names(criteria)) {
      for (o in outcome_vars) {
        adjmat <- get_adjacency(criterion, o, margin, bounds, diffs)
        plot_hasse(
          adjmat,
          main = paste0(criterion, "; ", outcome_labels[[o]], "; ", margin_label(margin)),
          horizontal_spacing = 2.5,
          vertical_spacing = 3,
          component_gap = 3
        )
      }
    }
    dev.off()
  }

  # Hurwicz over an alpha grid: one page per outcome and margin
  alpha_seq <- seq(0, 1, by = 0.1)
  n_col <- ceiling(sqrt(length(alpha_seq)))
  n_row <- ceiling(length(alpha_seq) / n_col)

  pdf(file.path(fig_dir, "Hurwicz_graphs.pdf"), width = 20, height = 15)
  for (o in outcome_vars) {
    for (margin in c(0, 2 / 12)) {
      par(mfrow = c(n_row, n_col), mar = c(1, 1, 3, 1), oma = c(0, 0, 3, 0))
      for (a in alpha_seq) {
        adjmat <- get_adjacency("Hurwicz", o, margin, bounds, diffs, alpha = a)
        plot_hasse(adjmat, main = paste0("alpha = ", a), vertex_size = 20, label_cex = 1.2)
      }
      mtext(paste0("Hurwicz; ", outcome_labels[[o]], "; ", margin_label(margin)),
            outer = TRUE, line = 1, cex = 1.5)
    }
  }
  dev.off()
}
