# Pareto front
#
# Mean OS / PFS against 5-year MDS/AML CIR for the two extreme scenarios, with
# 95% highest-density regions across the 1000 sampled parameter sets
# (probabilistic sensitivity analysis), and the same outcomes as differences
# from a reference strategy (PEC) within each parameter set.
# Input:  Data/simulation_outputs/<set>/NEC_data_uncertainty_bounds.rds
#         (from analysis/01_oscar_simulate_uncertainty_bounds.R)
# Output: Fig/<set>/Pareto_front_os_cloud.pdf, Pareto_front_pfs_cloud.pdf,
#         mean_outcomes_2scenarios.csv,
#         Pareto_difference_vs_PEC_os.pdf, Pareto_difference_vs_PEC_pfs.pdf,
#         differences_vs_PEC.csv, probability_best.csv

devtools::load_all()
library(data.table)
library(MASS)
library(ggplot2)
library(ggrepel)

# Parameter sets to plot; sets without simulation outputs are skipped
param_sets <- c("NETTER&RADIANT", "COMPETE")
input_root <- file.path("Data", "simulation_outputs")
fig_root <- "Fig"

mean_ci <- function(x, level = 0.95, na.rm = TRUE) {
  alpha <- 1 - level
  qs <- quantile(
    x,
    probs = c(alpha / 2, 1 - alpha / 2),
    na.rm = na.rm,
    names = FALSE
  )
  sprintf("%.2f (%.2f, %.2f)", mean(x, na.rm = na.rm), qs[1], qs[2])
}

# Kernel density estimation

#' Estimate a two-dimensional highest-density region
#'
#' @param x Numeric x values.
#' @param y Numeric y values.
#' @param prob Probability mass contained in the highest-density region.
#' @param n Integer grid resolution.
#' @param h Numeric vector of x- and y-axis KDE bandwidths.
#' @param pad Numeric bandwidth multiplier used to extend the KDE grid.
#'
#' @return A data.table containing the KDE grid and HDR indicators.
#'
#' @keywords internal
get_kde_hdr <- function(x,
                        y,
                        prob = 0.95,
                        n = 150L,
                        h = NULL,
                        pad = 3) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  
  if (data.table::uniqueN(x) < 2L ||
      data.table::uniqueN(y) < 2L) {
    stop("Both variables must contain at least two distinct finite values.")
  }
  
  if (is.null(h)) {
    h <- c(bandwidth.nrd(x), bandwidth.nrd(y))
  }
  
  if (length(h) != 2L || any(!is.finite(h)) || any(h <= 0)) {
    stop("`h` must contain two positive bandwidths.")
  }
  
  kd <- kde2d(
    x = x,
    y = y,
    n = n,
    h = h,
    lims = c(range(x) + c(-1, 1) * pad * h[1], range(y) + c(-1, 1) * pad * h[2])
  )
  
  z <- as.vector(kd$z)
  cell_mass <- z * diff(kd$x)[1L] * diff(kd$y)[1L]
  ord <- order(z, decreasing = TRUE)
  cumulative_mass <- cumsum(cell_mass[ord]) / sum(cell_mass)
  cutoff <- z[ord][which(cumulative_mass >= prob)[1L]]
  
  grid <- as.data.table(expand.grid(x = kd$x, y = kd$y))
  
  grid[, `:=`(density_ratio = z / cutoff,
              inside_hdr = z >= cutoff)]
  
  grid[]
}

# Strategy groups used for colour (position of PRRT) and line type (number of drugs)
prrt_position_levels <- c("PRRT first", "PRRT second", "PRRT third")
n_drugs_levels <- c("3-drug strategy", "2-drug strategy")

#' Add PRRT position and number of drugs for each strategy
#'
#' @param dt A data.table with a `tx` column of strategy names.
#'
#' @return `dt`, modified by reference, with factor columns `prrt_position`
#'   and `n_drugs`.
#'
#' @keywords internal
add_strategy_groups <- function(dt) {
  sequences <- treatment_sequences()[as.character(dt$tx)]
  position <- vapply(sequences, function(x) match("prrt", x), integer(1))
  n_drugs <- vapply(sequences, function(x) sum(!is.na(x)), integer(1))
  dt[, `:=`(
    prrt_position = factor(prrt_position_levels[position], levels = prrt_position_levels),
    n_drugs = factor(paste0(n_drugs, "-drug strategy"), levels = n_drugs_levels)
  )]
  dt[]
}

#' Add one highest-density region per scenario and treatment
#'
#' Region borders are coloured by PRRT position and drawn dashed for two-drug
#' strategies (mapped aesthetics, so they appear in the legend).
#'
#' @param p A ggplot object.
#' @param hdr_grid Output of get_kde_hdr() by scenario and tx, with the columns
#'   added by add_strategy_groups().
#' @param cloud_alpha Numeric transparency of filled regions.
#' @param style `"outline"` draws only the region's border; `"fill"` also
#'   fills it.
#' @param linewidth Numeric width of the border.
#' @param line_alpha Numeric transparency of the border.
#'
#' @return The ggplot object with the regions added.
#'
#' @keywords internal
add_hdr_clouds <- function(p, hdr_grid, cloud_alpha = 0.1,
                           style = c("outline", "fill"), linewidth = 0.6, line_alpha = 1) {
  style <- match.arg(style)
  # Contours must be added separately for each regular grid.
  groups <- unique(hdr_grid[, .(scenario, tx)])
  for (i in seq_len(nrow(groups))) {
    grid_i <- hdr_grid[scenario == groups$scenario[i] & tx == groups$tx[i]]
    if (!any(grid_i$inside_hdr, na.rm = TRUE)) {
      next
    }
    if (style == "fill") {
      p <- p +
        ggplot2::geom_contour_filled(
          data = grid_i,
          ggplot2::aes(x = x, y = y, z = density_ratio, fill = prrt_position),
          breaks = c(1, max(grid_i$density_ratio, na.rm = TRUE) + 1e-6),
          alpha = cloud_alpha,
          inherit.aes = FALSE,
          show.legend = FALSE
        )
    }
    # The region's border is the density_ratio = 1 contour
    p <- p +
      ggplot2::geom_contour(
        data = grid_i,
        ggplot2::aes(x = x, y = y, z = density_ratio, colour = prrt_position, linetype = n_drugs),
        breaks = 1,
        linewidth = if (style == "fill") 0.2 else linewidth,
        alpha = line_alpha,
        inherit.aes = FALSE
      )
  }
  p
}

#' Shared scales, legend and theme for the Pareto figures
#'
#' @param palette Named character vector of colours by PRRT position.
#'
#' @return A list of ggplot components.
#'
#' @keywords internal
pareto_scales_and_theme <- function(palette) {
  list(
    ggplot2::scale_colour_manual(values = palette, name = NULL, drop = FALSE),
    ggplot2::scale_fill_manual(values = palette, guide = "none", drop = FALSE),
    ggplot2::scale_linetype_manual(values = c("3-drug strategy" = "solid", "2-drug strategy" = "dashed"),
                                   name = NULL, drop = FALSE),
    # Legend keys opaque, whatever the border transparency
    ggplot2::guides(
      colour = ggplot2::guide_legend(order = 1, override.aes = list(alpha = 1, linewidth = 0.8)),
      linetype = ggplot2::guide_legend(order = 2, override.aes = list(alpha = 1, linewidth = 0.6,
                                                                      colour = "grey30"))
    ),
    ggplot2::theme_bw(),
    ggplot2::theme(
      legend.position = "bottom",
      legend.key.width = ggplot2::unit(1.5, "lines"),
      panel.grid = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 13)
    )
  )
}

#' Labels for every strategy, never dropped for overlap
#'
#' @param data Data with x_mean, y_mean, tx, prrt_position.
#' @param label_size Numeric label size.
#'
#' @return A ggrepel layer.
#'
#' @keywords internal
strategy_labels <- function(data, label_size) {
  ggrepel::geom_text_repel(
    data = data,
    ggplot2::aes(x = x_mean, y = y_mean, label = tx, colour = prrt_position),
    size = label_size,
    max.overlaps = Inf,
    min.segment.length = 0,
    box.padding = 0.4,
    seed = 1,
    show.legend = FALSE,
    inherit.aes = FALSE
  )
}

# Plot function

#' Plot treatment outcomes with joint uncertainty regions
#'
#' @description
#' Plots mean survival against 5-year MDS/AML cumulative incidence by
#' treatment sequence and scenario. Points show mean outcomes across
#' parameter sets, and regions show estimated joint highest-density
#' uncertainty regions, coloured by the position of PRRT.
#'
#' @param dat A data.table containing simulation results.
#' @param y_var Character string naming the survival outcome. Typically
#'   `"os_mean"` or `"pfs_mean"`.
#' @param palette Named character vector of colours by PRRT position.
#' @param prob Probability mass contained in each highest-density region.
#' @param h Numeric vector of KDE bandwidths for MDS/AML incidence and
#'   survival, respectively; `NULL` estimates them from the data.
#' @param n Integer KDE grid resolution.
#' @param x_limits,y_limits Optional numeric vectors of axis display limits.
#' @param point_size Numeric size of mean points.
#' @param label_size Numeric size of treatment labels.
#' @param cloud_alpha Numeric transparency of filled regions.
#' @param cloud_style `"outline"` or `"fill"`.
#' @param line_alpha Numeric transparency of region borders.
#'
#' @return A ggplot object.
#'
#' @keywords internal
plot_pareto_cloud <- function(dat,
                              y_var = c("os_mean", "pfs_mean"),
                              palette,
                              prob = 0.95,
                              h = NULL,
                              n = 150L,
                              x_limits = NULL,
                              y_limits = NULL,
                              point_size = 1.8,
                              label_size = 3,
                              cloud_alpha = 0.1,
                              cloud_style = c("outline", "fill"),
                              line_alpha = 1) {
  y_var <- match.arg(y_var)
  cloud_style <- match.arg(cloud_style)
  required <- c("scenario", "tx", "mds_aml_5y_cir", y_var)
  
  if (!all(required %in% names(dat))) {
    stop("Missing required columns: ", paste(setdiff(required, names(dat)), collapse = ", "))
  }
  
  mean_points <- add_strategy_groups(dat[, .(
    x_mean = mean(mds_aml_5y_cir, na.rm = TRUE),
    y_mean = mean(get(y_var), na.rm = TRUE)
  ), by = .(scenario, tx)])
  
  hdr_grid <- add_strategy_groups(dat[, get_kde_hdr(
    x = mds_aml_5y_cir,
    y = get(y_var),
    prob = prob,
    n = n,
    h = h
  ), by = .(scenario, tx)])
  
  p <- add_hdr_clouds(ggplot2::ggplot(), hdr_grid, cloud_alpha,
                      style = cloud_style, line_alpha = line_alpha)
  
  p +
    ggplot2::geom_point(
      data = mean_points,
      ggplot2::aes(x = x_mean, y = y_mean, colour = prrt_position),
      size = point_size,
      show.legend = FALSE,
      inherit.aes = FALSE
    ) +
    strategy_labels(mean_points, label_size) +
    ggplot2::facet_wrap(~ scenario) +
    ggplot2::scale_x_continuous(labels = scales::label_percent(accuracy = 1)) +
    ggplot2::labs(
      x = "5-year MDS/AML CIR, %",
      y = switch(y_var, os_mean = "Mean OS, years", pfs_mean = "Mean PFS, years")
    ) +
    pareto_scales_and_theme(palette) +
    ggplot2::theme(axis.title = ggplot2::element_text(size = 14)) +
    ggplot2::coord_cartesian(xlim = x_limits, ylim = y_limits)
}

#' Differences from a reference strategy within each parameter set
#'
#' @description
#' All strategies are simulated with the same parameter sets, so outcomes can
#' be compared set by set. Uncertainty shared by all strategies cancels out,
#' leaving the uncertainty in the comparison itself.
#'
#' @param dat A data.table of simulation results with `set_id`, `scenario`,
#'   `tx`, `os_mean`, `pfs_mean`, and `mds_aml_5y_cir`.
#' @param reference Character string naming the reference strategy.
#'
#' @return A data.table with one row per parameter set, scenario, and
#'   non-reference strategy, with differences in mean OS and PFS (months) and
#'   in 5-year MDS/AML CIR (percentage points).
#'
#' @keywords internal
compute_differences <- function(dat, reference = "PCE") {
  stopifnot(reference %in% dat$tx)
  ref <- dat[tx == reference, .(set_id, scenario,
                                os_ref = os_mean, pfs_ref = pfs_mean, cir_ref = mds_aml_5y_cir)]
  diffs <- ref[dat[tx != reference], on = .(set_id, scenario)]
  diffs[, `:=`(
    os_mean_diff_months = (os_mean - os_ref) * 12,
    pfs_mean_diff_months = (pfs_mean - pfs_ref) * 12,
    mds_aml_5y_cir_diff_pp = (mds_aml_5y_cir - cir_ref) * 100
  )]
  diffs[, tx := droplevels(tx)]
  diffs[]
}

#' Summarise differences from the reference strategy
#'
#' @param diffs Output of compute_differences().
#'
#' @return A data.table by strategy and scenario with the mean (95% interval)
#'   of each difference and the share of parameter sets in each quadrant.
#'
#' @keywords internal
summarise_differences <- function(diffs) {
  diffs[, .(
    os_mean_diff_months = mean_ci(os_mean_diff_months),
    pfs_mean_diff_months = mean_ci(pfs_mean_diff_months),
    mds_aml_5y_cir_diff_pp = mean_ci(mds_aml_5y_cir_diff_pp),
    p_longer_os = mean(os_mean_diff_months > 0),
    p_lower_cir = mean(mds_aml_5y_cir_diff_pp < 0),
    # Quadrants of the OS figure
    p_longer_os_lower_cir = mean(os_mean_diff_months > 0 & mds_aml_5y_cir_diff_pp < 0),
    p_longer_os_higher_cir = mean(os_mean_diff_months > 0 & mds_aml_5y_cir_diff_pp >= 0),
    p_shorter_os_lower_cir = mean(os_mean_diff_months <= 0 & mds_aml_5y_cir_diff_pp < 0),
    p_shorter_os_higher_cir = mean(os_mean_diff_months <= 0 & mds_aml_5y_cir_diff_pp >= 0)
  ), by = .(tx, scenario)][order(scenario, tx)]
}

#' Share of parameter sets in which each strategy has the longest mean survival
#'
#' @param dat A data.table of simulation results.
#'
#' @return A data.table by strategy and scenario.
#'
#' @keywords internal
probability_best <- function(dat) {
  best <- function(var) {
    winners <- dat[, .SD[which.max(get(var))], by = .(set_id, scenario)][, .(set_id, scenario, tx)]
    winners[, .(p = .N / uniqueN(dat$set_id)), by = .(scenario, tx)]
  }
  grid <- CJ(scenario = unique(dat$scenario), tx = unique(dat$tx))
  out <- merge(grid, best("os_mean")[, .(scenario, tx, p_best_os = p)], all.x = TRUE)
  out <- merge(out, best("pfs_mean")[, .(scenario, tx, p_best_pfs = p)], all.x = TRUE)
  out[is.na(p_best_os), p_best_os := 0]
  out[is.na(p_best_pfs), p_best_pfs := 0]
  out[order(scenario, tx)]
}

#' Plot differences from the reference strategy with joint uncertainty regions
#'
#' @param diffs Output of compute_differences().
#' @param y_var `"os_mean_diff_months"` or `"pfs_mean_diff_months"`.
#' @param reference Character string naming the reference strategy.
#' @param palette Named character vector of colours by PRRT position.
#' @param prob Probability mass contained in each highest-density region.
#' @param n Integer KDE grid resolution.
#' @param x_limits,y_limits Optional numeric vectors of axis display limits.
#' @param point_size Numeric size of mean points.
#' @param label_size Numeric size of treatment labels.
#' @param cloud_alpha Numeric transparency of filled regions.
#' @param cloud_style `"outline"` or `"fill"`.
#' @param line_alpha Numeric transparency of region borders.
#'
#' @return A ggplot object.
#'
#' @keywords internal
plot_difference_cloud <- function(diffs,
                                  y_var = c("os_mean_diff_months", "pfs_mean_diff_months"),
                                  reference = "PCE",
                                  palette,
                                  prob = 0.95,
                                  n = 150L,
                                  x_limits = NULL,
                                  y_limits = NULL,
                                  point_size = 1.8,
                                  label_size = 3,
                                  cloud_alpha = 0.1,
                                  cloud_style = c("outline", "fill"),
                                  line_alpha = 1) {
  y_var <- match.arg(y_var)
  cloud_style <- match.arg(cloud_style)
  outcome_label <- switch(y_var, os_mean_diff_months = "OS", pfs_mean_diff_months = "PFS")
  
  mean_points <- add_strategy_groups(diffs[, .(
    x_mean = mean(mds_aml_5y_cir_diff_pp),
    y_mean = mean(get(y_var))
  ), by = .(scenario, tx)])
  # The reference strategy sits at the origin
  reference_point <- add_strategy_groups(
    unique(diffs[, .(scenario)])[, `:=`(tx = reference, x_mean = 0, y_mean = 0)]
  )
  
  hdr_grid <- add_strategy_groups(diffs[, get_kde_hdr(
    x = mds_aml_5y_cir_diff_pp,
    y = get(y_var),
    prob = prob,
    n = n
  ), by = .(scenario, tx)])
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3)
  p <- add_hdr_clouds(p, hdr_grid, cloud_alpha, style = cloud_style, line_alpha = line_alpha)
  
  p +
    ggplot2::geom_point(
      data = mean_points,
      ggplot2::aes(x = x_mean, y = y_mean, colour = prrt_position),
      size = point_size,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = reference_point,
      ggplot2::aes(x = x_mean, y = y_mean, colour = prrt_position),
      size = point_size * 1.5,
      show.legend = FALSE
    ) +
    strategy_labels(rbind(mean_points, reference_point, use.names = TRUE), label_size) +
    ggplot2::facet_wrap(~ scenario) +
    ggplot2::scale_y_continuous(breaks = scales::breaks_width(5)) +
    ggplot2::labs(
      x = paste0("Difference in 5-year MDS/AML CIR vs ", reference, ", percentage points"),
      y = paste0("Difference in mean ", outcome_label, " vs ", reference, ", months")
    ) +
    pareto_scales_and_theme(palette) +
    ggplot2::theme(axis.title = ggplot2::element_text(size = 12)) +
    ggplot2::coord_cartesian(xlim = x_limits, ylim = y_limits)
}

# One colour per position of PRRT in the sequence; two-drug strategies are
# drawn with dashed borders (see the legend).
cb_palette <- c(
  "PRRT first" = "#2171B5",   # blue:   PEC, PCE, PE
  "PRRT second" = "#D94801",  # orange: EPC, CPE, EP
  "PRRT third" = "#6A51A3"    # purple: ECP, CEP
)
# "outline" (region borders only) or "fill"
cloud_style <- "outline"
# Transparency of the region borders (1 = opaque)
line_alpha <- 0.3

for (param_set in param_sets) {
  input_file <- file.path(input_root, param_set, "NEC_data_uncertainty_bounds.rds")
  if (!file.exists(input_file)) {
    message("Skipping ", param_set, ": ", input_file, " not found")
    next
  }
  message("=== Parameter set: ", param_set)
  fig_dir <- file.path(fig_root, param_set)
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Mean OS, PFS, MDS/AML CIR and 95% intervals across the parameter sets
  dat <- readRDS(input_file)
  dat[, scenario := factor(
    fifelse(
      discount1 == 1 & discount2 == 1 & aml_mds_multiplier == 1,
      "No Tx Discount, Low MDS/AML",
      "High Tx Discount, High MDS/AML"
    ),
    levels = c("No Tx Discount, Low MDS/AML", "High Tx Discount, High MDS/AML")
  )]
  dat[, tx := factor(tx, levels = names(treatment_sequences()))]
  
  mean_outcomes_2scenarios <- dat[, .(
    os_mean = mean_ci(os_mean),
    pfs_mean = mean_ci(pfs_mean),
    mds_aml_5y_cir_pct = mean_ci(mds_aml_5y_cir * 100)
  ), by = .(tx, scenario)]
  setorder(mean_outcomes_2scenarios, scenario, tx)
  write.csv(mean_outcomes_2scenarios, file.path(fig_dir, "mean_outcomes_2scenarios.csv"), row.names = FALSE)
  
  for (y_var in c("os_mean", "pfs_mean")) {
    # Axes fitted by ggplot to this figure's regions and points
    p <- plot_pareto_cloud(
      dat = dat,
      y_var = y_var,
      palette = cb_palette,
      h = NULL, # bandwidths from the data (normal reference rule), per cloud
      label_size = 4,
      cloud_style = cloud_style,
      line_alpha = line_alpha
    )
    ggplot2::ggsave(
      filename = file.path(fig_dir, paste0("Pareto_front_", sub("_mean$", "", y_var), "_cloud.pdf")),
      plot = p,
      width = 7,
      height = 5
    )
  }
  
  # Differences from the reference strategy within each parameter set
  reference <- "PCE"
  diffs <- compute_differences(dat, reference = reference)
  write.csv(summarise_differences(diffs),
            file.path(fig_dir, paste0("differences_vs_", reference, ".csv")), row.names = FALSE)
  write.csv(probability_best(dat), file.path(fig_dir, "probability_best.csv"), row.names = FALSE)
  
  for (y_var in c("os_mean_diff_months", "pfs_mean_diff_months")) {
    # Axes fitted by ggplot to this figure's regions, points and the origin
    p <- plot_difference_cloud(
      diffs = diffs,
      y_var = y_var,
      reference = reference,
      palette = cb_palette,
      label_size = 4,
      cloud_style = cloud_style,
      line_alpha = line_alpha
    )
    ggplot2::ggsave(
      filename = file.path(fig_dir, paste0("Pareto_difference_vs_", reference, "_",
                                           sub("_mean_diff_months$", "", y_var), ".pdf")),
      plot = p,
      width = 8,
      height = 5
    )
  }
}
