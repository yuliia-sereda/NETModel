#' Outcomes reported in the clinician app
#'
#' @description
#' One row per outcome produced by [summarise_population_app()], with its
#' label, units, the scaling from simulated values to reported values and
#' differences, and whether higher values are better.
#'
#' @return A `data.table` with columns `outcome`, `label`, `unit`, `scale`,
#'   `diff_unit`, `diff_scale`, and `higher_is_better`.
#'
#' @export
app_outcome_info <- function() {
  survival <- function(prefix, name) {
    data.table(
      outcome = paste0(prefix, c("_mean", "_median", "_1y", "_3y", "_5y")),
      label = c(paste("Mean", name), paste("Median", name), paste(name, "at", c(1, 3, 5), "years")),
      unit = c("years", "years", "%", "%", "%"),
      scale = c(1, 1, 100, 100, 100),
      diff_unit = c("months", "months", "percentage points", "percentage points", "percentage points"),
      diff_scale = c(12, 12, 100, 100, 100),
      higher_is_better = TRUE
    )
  }
  rbind(
    survival("os", "OS"),
    survival("pfs", "PFS"),
    data.table(
      outcome = paste0("mds_aml_", c(1, 3, 5), "y_cir"),
      label = paste("MDS/AML cumulative incidence at", c(1, 3, 5), "years"),
      unit = "%",
      scale = 100,
      diff_unit = "percentage points",
      diff_scale = 100,
      higher_is_better = FALSE
    )
  )
}

#' Summarise app outcomes across parameter sets
#'
#' @description
#' For each design point, treatment sequence, and outcome: the mean and the
#' 95% uncertainty interval (2.5th and 97.5th percentiles) across parameter
#' sets, in reporting units (years or %).
#'
#' @param sets A `data.table` of set-level results from [sim_parameter_sets()]
#'   with `summarise = summarise_population_app`.
#' @param design Optional `data.table` of design points with a `design_id`
#'   column; its columns are joined to the summary.
#'
#' @return A long `data.table` with one row per design point, treatment
#'   sequence, and outcome.
#'
#' @export
summarise_app_outcomes <- function(sets, design = NULL) {
  # Column names used in data.table expressions (declared for R CMD check)
  tx <- design_id <- outcome <- label <- unit <- value <- NULL
  
  info <- app_outcome_info()
  long <- melt(sets, id.vars = c("set_id", "design_id", "tx"), measure.vars = info$outcome,
               variable.name = "outcome", value.name = "value", variable.factor = FALSE)
  long <- info[, list(outcome, label, unit, scale)][long, on = "outcome"]
  out <- long[, list(
    mean = mean(value * scale),
    lower95 = stats::quantile(value * scale, 0.025, names = FALSE),
    upper95 = stats::quantile(value * scale, 0.975, names = FALSE),
    n_sets = .N
  ), by = list(design_id, tx, outcome, label, unit)]
  if (!is.null(design)) {
    out <- design[out, on = "design_id"]
  }
  out[]
}

#' Pairwise differences between treatment sequences across parameter sets
#'
#' @description
#' For each design point, pair of treatment sequences (tx1, tx2), and outcome:
#' the difference tx1 - tx2 computed within each parameter set, summarised by
#' its mean, 95% uncertainty interval, and the share of parameter sets in
#' which tx1 is better (higher survival, or lower MDS/AML incidence).
#' Differences are in months for mean and median survival and in percentage
#' points for probabilities.
#'
#' @inheritParams summarise_app_outcomes
#'
#' @return A long `data.table` with one row per design point, pair, and
#'   outcome.
#'
#' @export
summarise_app_differences <- function(sets, design = NULL) {
  # Column names used in data.table expressions (declared for R CMD check)
  design_id <- outcome <- label <- diff_unit <- diff_scale <- higher_is_better <- NULL
  
  info <- app_outcome_info()
  long <- melt(sets, id.vars = c("set_id", "design_id", "tx"), measure.vars = info$outcome,
               variable.name = "outcome", value.name = "value", variable.factor = FALSE)
  wide <- dcast(long, set_id + design_id + outcome ~ tx, value.var = "value")
  wide <- info[, list(outcome, diff_unit, diff_scale, higher_is_better)][wide, on = "outcome"]

  tx_names <- intersect(names(treatment_sequences()), unique(as.character(sets$tx)))
  pairs <- utils::combn(tx_names, 2, simplify = FALSE)
  out <- rbindlist(lapply(pairs, function(pair) {
    d <- (wide[[pair[1]]] - wide[[pair[2]]]) * wide$diff_scale
    better <- data.table::fifelse(wide$higher_is_better, d > 0, d < 0)
    data.table(design_id = wide$design_id, outcome = wide$outcome, diff_unit = wide$diff_unit,
               d = d, better = better)[, list(
      mean_diff = mean(d),
      lower95 = stats::quantile(d, 0.025, names = FALSE),
      upper95 = stats::quantile(d, 0.975, names = FALSE),
      p_tx1_better = mean(better),
      n_sets = .N
    ), by = list(design_id, outcome, diff_unit)][, `:=`(tx1 = pair[1], tx2 = pair[2])]
  }))
  out <- info[, list(outcome, label)][out, on = "outcome"]
  setcolorder(out, c("design_id", "tx1", "tx2", "outcome", "label", "diff_unit"))
  if (!is.null(design)) {
    out <- design[out, on = "design_id"]
  }
  out[]
}
