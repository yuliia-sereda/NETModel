# Helpers for the clinician app. The app only reads precomputed simulation
# outputs from data/; it never runs the model.

# Shiny sources this file before app.R, so load what the helpers use here
library(data.table)

# ------------------------------------------------------------------------------
# Strategies, ordered by the position of PRRT in the sequence
# ------------------------------------------------------------------------------
position_levels <- c("PRRT first", "PRRT second", "PRRT third")
drugs_levels <- c("3-drug strategy", "2-drug strategy")

strategies <- data.table(
  tx = c(
    "PE", 
    "EP", 
    "PEC", 
    "PCE", 
    "EPC", 
    "ECP", 
    "CPE", 
    "CEP"),
  name = c(
    "PRRT → everolimus",
    "everolimus → PRRT",
    "PRRT → everolimus → cabozantinib",
    "PRRT → cabozantinib → everolimus",
    "everolimus → PRRT → cabozantinib",
    "everolimus → cabozantinib → PRRT",
    "cabozantinib → PRRT → everolimus",
    "cabozantinib → everolimus → PRRT"
  ),
  prrt_position = factor(position_levels[c(1, 2, 1, 1, 2, 3, 2, 3)], levels = position_levels),
  n_drugs = factor(c(
    rep("2-drug strategy", 2),
    rep("3-drug strategy", 6)
    ), levels = drugs_levels)
)
strategy_order <- strategies$tx
strategy_choices <- setNames(strategies$tx, paste0(strategies$tx, ": ", strategies$name))

position_palette <- c("PRRT first" = "#2171B5", "PRRT second" = "#D94801", "PRRT third" = "#6A51A3")
drugs_linetypes <- c("3-drug strategy" = "solid", "2-drug strategy" = "dashed")

add_strategy_info <- function(dt) {
  strategies[, .(tx, name, prrt_position, n_drugs)][dt, on = "tx"]
}

# ------------------------------------------------------------------------------
# Outcomes
# ------------------------------------------------------------------------------
outcome_choices <- function(info) {
  split_info <- split(info, factor(info$group, levels = unique(info$group)))
  lapply(split_info, function(x) setNames(x$outcome, x$label))
}

# Years with 2 decimals; months, % and percentage points with 1. `unit` may be a
# single unit or one per value.
format_value <- function(x, unit) {
  fmt <- ifelse(rep_len(unit, length(x)) == "years", "%.2f", "%.1f")
  sprintf(fmt, x)
}

format_estimate <- function(est, lower, upper, unit) {
  paste0(format_value(est, unit), " (", format_value(lower, unit), "–",
         format_value(upper, unit), ")")
}

unit_suffix <- function(unit) {
  switch(unit, years = "years", "%" = "%", months = "months", "percentage points" = "percentage points")
}

# ------------------------------------------------------------------------------
# 95% highest-density region of a cloud of points, as contour paths
# ------------------------------------------------------------------------------
hdr_contour <- function(x, y, prob = 0.95, n = 100, pad = 3) {
  h <- c(MASS::bandwidth.nrd(x), MASS::bandwidth.nrd(y))
  h[!is.finite(h) | h <= 0] <- 1e-6
  kd <- MASS::kde2d(x, y, n = n, h = h,
                    lims = c(range(x) + c(-1, 1) * pad * h[1], range(y) + c(-1, 1) * pad * h[2]))
  cell_mass <- kd$z * diff(kd$x[1:2]) * diff(kd$y[1:2])
  ord <- order(kd$z, decreasing = TRUE)
  cutoff <- kd$z[ord][which(cumsum(cell_mass[ord]) / sum(cell_mass) >= prob)[1]]
  lines <- grDevices::contourLines(kd$x, kd$y, kd$z, levels = cutoff)
  rbindlist(lapply(seq_along(lines), function(i) data.table(x = lines[[i]]$x, y = lines[[i]]$y, piece = i)))
}

# ------------------------------------------------------------------------------
# Plots (all strategies, coloured by PRRT position, dashed for 2-drug strategies)
# ------------------------------------------------------------------------------
# Scales for the aesthetics a plot maps (ggplot2 warns about unused manual scales)
strategy_scales <- function(aesthetics = c("colour", "linetype")) {
  scales <- list(
    colour = ggplot2::scale_colour_manual(values = position_palette, name = NULL, drop = FALSE),
    fill = ggplot2::scale_fill_manual(values = position_palette, name = NULL, drop = FALSE),
    linetype = ggplot2::scale_linetype_manual(values = drugs_linetypes, name = NULL, drop = FALSE)
  )
  c(
    unname(scales[aesthetics]),
    list(ggplot2::theme_bw(base_size = 14),
    ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank(),
                   legend.key.width = ggplot2::unit(1.6, "lines")))
  )
}

strategy_labels <- function(data, size = 4.5) {
  ggrepel::geom_text_repel(
    data = data, ggplot2::aes(x = x, y = y, label = tx, colour = prrt_position),
    size = size, max.overlaps = Inf, min.segment.length = 0, box.padding = 0.4,
    seed = 1, show.legend = FALSE
  )
}

# `cloud` has one row per parameter set and strategy, with columns tx, x, y.
plot_cloud <- function(cloud, x_label, y_label, reference = NULL) {
  regions <- cloud[, hdr_contour(x, y), by = tx]
  regions[, group := paste(tx, piece)]
  regions <- add_strategy_info(regions)
  means <- add_strategy_info(cloud[, .(x = mean(x), y = mean(y)), by = tx])
  p <- ggplot2::ggplot()
  if (!is.null(reference)) {
    p <- p +
      ggplot2::geom_hline(yintercept = 0, colour = "grey60") +
      ggplot2::geom_vline(xintercept = 0, colour = "grey60")
    ref_point <- add_strategy_info(data.table(tx = reference, x = 0, y = 0))
    means <- rbind(means, ref_point, use.names = TRUE)
  }
  p +
    ggplot2::geom_path(data = regions, ggplot2::aes(x, y, group = group, colour = prrt_position,
                                                    linetype = n_drugs), linewidth = 0.7, alpha = 0.55) +
    ggplot2::geom_point(data = means, ggplot2::aes(x, y, colour = prrt_position), size = 2.6,
                        show.legend = FALSE) +
    strategy_labels(means) +
    ggplot2::labs(x = x_label, y = y_label) +
    strategy_scales() +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(alpha = 1, linewidth = 1)),
                    linetype = ggplot2::guide_legend(override.aes = list(alpha = 1, colour = "grey30")))
}

# Point estimate and 95% range for each strategy
# `show_point = FALSE` draws only the range (lower95 to upper95), e.g. for the
# worst-to-best case range of a strategy, where a midpoint has no meaning.
plot_intervals <- function(dt, x_label, show_point = TRUE) {
  dt <- add_strategy_info(copy(dt))
  dt[, tx := factor(tx, levels = rev(strategy_order))]
  p <- ggplot2::ggplot(dt, ggplot2::aes(y = tx, xmin = lower95, xmax = upper95,
                                        colour = prrt_position, linetype = n_drugs)) +
    ggplot2::geom_linerange(linewidth = 1.2)
  if (show_point) {
    p <- p + ggplot2::geom_point(ggplot2::aes(x = mean), size = 3)
  }
  p +
    ggplot2::labs(x = x_label, y = NULL) +
    strategy_scales()
}

# ------------------------------------------------------------------------------
# Comparisons and choice functions
# ------------------------------------------------------------------------------

# Differences from a reference strategy within each parameter set, for the
# uncertainty regions: x = MDS/AML at 5 years (percentage points),
# y = mean OS or PFS (months).
difference_cloud <- function(cloud, reference, y_var) {
  ref <- cloud[tx == reference, .(set_id, ref_y = get(y_var), ref_x = mds_aml_5y_cir)]
  ref[cloud[tx != reference], on = "set_id"][
    , .(tx = as.character(tx), x = (mds_aml_5y_cir - ref_x) * 100, y = (get(y_var) - ref_y) * 12)]
}

# Optimal strategies over a set of design points (combinations of assumptions)
#
# `means` has one row per design point and strategy with the average outcome
# (column `mean`). For each strategy:
#   worst, best   its least and most favourable average across the design points
#   maximality    not uniformly worse: no other strategy is better by more than
#                 `margin` at every design point
#   beaten_by     the strategies that are better at every design point
#   pessimism     best worst case (within `margin` of the highest worst case)
#   optimism      best best case (within `margin` of the highest best case)
# `margin` is in the units of `mean`; `higher_is_better` is FALSE for MDS/AML.
optimal_strategies <- function(means, higher_is_better = TRUE, margin = 0) {
  means <- copy(means)
  means[, score := if (higher_is_better) mean else -mean]
  wide <- dcast(means, design_id ~ tx, value.var = "score")
  txs <- intersect(strategy_order, names(wide))
  beaten_by <- vapply(txs, function(s) {
    others <- setdiff(txs, s)
    paste(others[vapply(others, function(t) min(wide[[t]] - wide[[s]]) > margin, logical(1))],
          collapse = ", ")
  }, character(1))
  ranges <- means[, .(worst = if (higher_is_better) min(mean) else max(mean),
                      best = if (higher_is_better) max(mean) else min(mean),
                      worst_score = min(score), best_score = max(score)), by = tx]
  ranges[, `:=`(
    maximality = beaten_by[tx] == "",
    beaten_by = beaten_by[tx],
    pessimism = worst_score >= max(worst_score) - margin,
    optimism = best_score >= max(best_score) - margin
  )]
  ranges[, tx := factor(tx, levels = strategy_order)]
  setorder(ranges, tx)
  ranges[, c("worst_score", "best_score") := NULL]
  ranges[]
}
