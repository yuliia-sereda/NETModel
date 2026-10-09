######################################
# Calibrate lambda_46 & lambda_dot5 #
#####################################

# This code calibrates lambda_46 and lambda_dot5 to match 
# the median OS in the NETTER-1 trial [48.0 months (95% CI 37.4-55.2)] and 
# 5-year CIR of MDS/AML in the Washington University in St. Louis data [0.051 (95% CI: 0.018-0.084)]
# (PE sequence: PRRT, then everolimus).

devtools::load_all()
library(data.table)
library(ggplot2)
library(lhs)
library(mgcv)

# Make calibration design
make_calibration_design <- function(
    n_design = 1000,
    seed = 20260728
) {
  set.seed(seed)
  
  valid_designs <- list()
  n_valid <- 0L
  
  while (n_valid < n_design) {
    
    # Generate more points than needed because some will be rejected
    n_batch <- max(100, 2 * (n_design - n_valid))
    
    X <- lhs::maximinLHS(n = n_batch, k = 4)
    
    batch <- data.table::data.table(

      # Median salvage-to-death time (1 week to 24 months, expressed in years)
      median_salvage_death_years =
        exp(log(1 / 52) +
              X[, 1] * (log(2) - log(1 / 52))),
      
      # Relative SE of the estimated median
      sdlog_salvage_death =
        0.05 + X[, 2] * (0.30 - 0.05),
      
      # Baseline MDS/AML hazard (0.001-0.010; 0.010 gives a 5-year CIR of ~13%)
      lambda_mdsaml =
        exp(log(0.001) +
              X[, 3] * (log(0.010) - log(0.001))),
      
      # SE of log(MDS/AML incidence rate)
      sdlog_lambda_mdsaml =
        0.05 + X[, 4] * (0.30 - 0.05)
    )
    
    # Convert relative SE to SE of the median survival time
    batch[, se_salvage_death_years :=
            median_salvage_death_years * sdlog_salvage_death]
    
    # Internal quality-control filters
    lower_tail_survival <- stats::qlnorm(
      0.001,
      meanlog = log(batch$median_salvage_death_years),
      sdlog = batch$sdlog_salvage_death
    )
    
    upper_tail_mdsaml <- stats::qlnorm(
      0.999,
      meanlog = log(batch$lambda_mdsaml),
      sdlog = batch$sdlog_lambda_mdsaml
    )
    
    batch <- batch[
      lower_tail_survival >= 1 / 365 &   # ≥ approximately 1 day
        upper_tail_mdsaml <= 0.10
    ]
    
    if (nrow(batch) > 0) {
      valid_designs[[length(valid_designs) + 1]] <- batch
      n_valid <- n_valid + nrow(batch)
    }
  }
  
  design <- data.table::rbindlist(valid_designs)[1:n_design]
  
  design[, run_id := .I]
  
  data.table::setcolorder(
    design,
    c(
      "run_id",
      "median_salvage_death_years",
      "se_salvage_death_years",
      "lambda_mdsaml",
      "sdlog_lambda_mdsaml"
    )
  )
  
  design[]
}

# Calibration function
calibrate_lambda_46_and_dot5 <- function(
    calibration_grid,
    param_input,
    seed,
    population,
    treatments,
    n_cores = 1L
) {
  D <- data.table::data.table(discount1 = 1, discount2 = 1, aml_mds_multiplier = 1)

  # lambda_46 and lambda_dot5 come from the calibration design only: any rows
  # for them in param_input (e.g. blank or earlier estimates) are ignored.
  calibrated_param_names <- c("lambda_46", "lambda_dot5")
  param_input_fixed <- param_input[!param %in% calibrated_param_names]

  run_calibration_point <- function(i) {
    param_input_i <- rbind(
      param_input_fixed,
      data.table::data.table(
        param = calibrated_param_names,
        point_stat = c("median", "ir"),
        est = c(calibration_grid[i, median_salvage_death_years], calibration_grid[i, lambda_mdsaml]),
        se = c(calibration_grid[i, se_salvage_death_years], calibration_grid[i, sdlog_lambda_mdsaml])
      ),
      fill = TRUE
    )

    NEC_data <- sim_data(
      param_input = param_input_i,
      design = D,
      population = population,
      treatments = treatments,
      seed = seed,
      show_progress = FALSE
    )
    
    NEC_data[, .(os_median, mds_aml_5y_cir)]
  }
  
  n_runs <- nrow(calibration_grid)
  
  if (n_cores > 1) {
    results <- parallel::mclapply(
      seq_len(n_runs),
      run_calibration_point,
      mc.cores = n_cores,
      mc.preschedule = FALSE
    )
    # A worker that errors returns a try-error; one that is killed (e.g. out of
    # memory) returns NULL. Treat both as failures rather than dropping rows.
    failed <- vapply(results, function(x) is.null(x) || inherits(x, "try-error"), logical(1))
    if (any(failed)) {
      stop("Calibration point(s) ", paste(which(failed), collapse = ", "), " failed: ",
           results[[which(failed)[1]]])
    }
  } else {
    pb <- utils::txtProgressBar(min = 0, max = n_runs, style = 3)
    on.exit(close(pb), add = TRUE)
    results <- vector("list", n_runs)
    for (i in seq_len(n_runs)) {
      results[[i]] <- run_calibration_point(i)
      utils::setTxtProgressBar(pb, i)
    }
  }
  
  cbind(calibration_grid, data.table::rbindlist(results))
}

# Calibration targets: NETTER-1 median OS and WashU 5-year CIR of MDS/AML
targets_dt <- data.table(
  param = c("os_median", "mds_aml_5y_cir"),
  target = c(48.0/12, 0.051),
  lower = c(37.4/12, 0.018),
  upper = c(55.2/12, 0.084)
)

os_median_target <- targets_dt[param == "os_median", target]
mds_aml_5y_cir_target <- targets_dt[param == "mds_aml_5y_cir", target]

# Emulator-based calibration
# Fit smooth emulators of both outputs over the calibration design and solve
# for the parameter values at which they hit the targets. Unlike picking the
# top-ranked design point, the emulators average over Monte Carlo noise.
# Median OS is driven mainly by lambda_46 and the 5-year CIR by lambda_dot5.
# The sdlog parameters have little effect on either output and cannot be
# calibrated, so they are held fixed at the values below (the same for every
# parameter set). They are also used as the se in the calibrated output:
# lambda_46 se = est * sdlog_salvage_death_fixed (se on the time scale), and
# lambda_dot5 se = sdlog_lambda_mdsaml_fixed (SD of the log rate).
sdlog_salvage_death_fixed <- 0.0685    # relative SE of the median salvage-to-death time
sdlog_lambda_mdsaml_fixed <- 0.278413  # SD of log(MDS/AML incidence rate)

fit_calibration_emulators <- function(calibration_results, set_name) {
  emulator_formula_rhs <- "~ s(log(median_salvage_death_years)) + s(log(lambda_mdsaml)) +
    s(sdlog_salvage_death) + s(sdlog_lambda_mdsaml)"
  emulator_os <- gam(as.formula(paste("os_median", emulator_formula_rhs)), data = calibration_results)
  emulator_cir <- gam(as.formula(paste("mds_aml_5y_cir", emulator_formula_rhs)), data = calibration_results)
  
  for (emulator in list(emulator_os, emulator_cir)) {
    message(sprintf(
      "Emulator for %s: deviance explained %.1f%%, residual SD %.4f",
      all.vars(formula(emulator))[1],
      100 * summary(emulator)$dev.expl,
      sqrt(emulator$sig2)
    ))
  }
  
  predict_emulators <- function(median_salvage_death_years, lambda_mdsaml) {
    newdata <- data.table(
      median_salvage_death_years = median_salvage_death_years,
      lambda_mdsaml = lambda_mdsaml,
      sdlog_salvage_death = sdlog_salvage_death_fixed,
      sdlog_lambda_mdsaml = sdlog_lambda_mdsaml_fixed
    )
    c(
      os_median = unname(predict(emulator_os, newdata)),
      mds_aml_5y_cir = unname(predict(emulator_cir, newdata))
    )
  }
  
  # Solve both emulators jointly (on the log scale, within the design range)
  solve_emulators <- function(os_median_goal, mds_aml_5y_cir_goal) {
    loss <- function(x) {
      pred <- predict_emulators(exp(x[1]), exp(x[2]))
      ((pred[["os_median"]] - os_median_goal) / os_median_goal)^2 +
        ((pred[["mds_aml_5y_cir"]] - mds_aml_5y_cir_goal) / mds_aml_5y_cir_goal)^2
    }
    limits <- log(rbind(
      range(calibration_results$median_salvage_death_years),
      range(calibration_results$lambda_mdsaml)
    ))
    fit <- optim(
      par = rowMeans(limits),
      fn = loss,
      method = "L-BFGS-B",
      lower = limits[, 1],
      upper = limits[, 2]
    )
    pred <- predict_emulators(exp(fit$par[1]), exp(fit$par[2]))
    data.table(
      os_median_goal = os_median_goal,
      mds_aml_5y_cir_goal = mds_aml_5y_cir_goal,
      median_salvage_death_years = exp(fit$par[1]),
      lambda_mdsaml = exp(fit$par[2]),
      os_median_emulated = pred[["os_median"]],
      mds_aml_5y_cir_emulated = pred[["mds_aml_5y_cir"]],
      relative_distance = sqrt(fit$value),
      at_design_bound = any(abs(fit$par - limits) < 1e-6)
    )
  }
  
  # Point estimate, plus the parameter values at each end of the target CIs
  # (varying one target at a time, holding the other at its point estimate)
  os_median_ci <- targets_dt[param == "os_median", c(lower, upper)]
  mds_aml_5y_cir_ci <- targets_dt[param == "mds_aml_5y_cir", c(lower, upper)]
  calibration_solutions <- rbindlist(list(
    target = solve_emulators(os_median_target, mds_aml_5y_cir_target),
    os_median_lower = solve_emulators(os_median_ci[1], mds_aml_5y_cir_target),
    os_median_upper = solve_emulators(os_median_ci[2], mds_aml_5y_cir_target),
    mds_aml_5y_cir_lower = solve_emulators(os_median_target, mds_aml_5y_cir_ci[1]),
    mds_aml_5y_cir_upper = solve_emulators(os_median_target, mds_aml_5y_cir_ci[2])
  ), idcol = "scenario")
  
  if (any(calibration_solutions$at_design_bound)) {
    warning(set_name, ": emulator solution(s) on the edge of the calibration design (target not reachable ",
            "within the design range): ",
            paste(calibration_solutions[at_design_bound == TRUE, scenario], collapse = ", "),
            ". Widen the ranges in make_calibration_design() if these are needed.")
  }
  
  # Values for the input_params file
  calibrated <- calibration_solutions[scenario == "target"]
  calibrated_params <- data.table(
    param = c("lambda_46", "lambda_dot5"),
    est = c(calibrated$median_salvage_death_years, calibrated$lambda_mdsaml),
    se = c(calibrated$median_salvage_death_years * sdlog_salvage_death_fixed, sdlog_lambda_mdsaml_fixed)
  )
  
  list(solutions = calibration_solutions, calibrated_params = calibrated_params)
}

plot_calibration <- function(calibration_results, calibration_solutions, set_name) {
  y_max <- max(
    ceiling(calibration_results$os_median),
    targets_dt[param == "os_median", upper],
    na.rm = TRUE
  )
  y_min <- min(
    floor(calibration_results$os_median),
    targets_dt[param == "os_median", lower],
    na.rm = TRUE
  )
  
  ggplot(data = calibration_results, 
         aes(y = os_median, x = mds_aml_5y_cir)
         ) +
    geom_point(size = 0.5) +
    geom_vline(xintercept = targets_dt[param == "mds_aml_5y_cir", c(lower, upper)],
               linetype = "dashed", color = "darkgrey") +
    geom_vline(xintercept = mds_aml_5y_cir_target,
               linetype = "solid", color = "darkgrey") +
    geom_hline(yintercept = targets_dt[param == "os_median", c(lower, upper)],
               linetype = "dashed", color = "darkgrey") +
    geom_hline(yintercept = os_median_target,
               linetype = "solid", color = "darkgrey") +
    geom_point(data = calibration_solutions[scenario == "target"],
               aes(x = mds_aml_5y_cir_emulated, y = os_median_emulated),
               color = "red", size = 2) +
    scale_y_continuous(limits = c(y_min, y_max),
                       breaks = seq(y_min, y_max, by = 1)
                       ) + 
    labs(y = "Median OS",
         x = "MDS/AML 5-year CIR",
         title = paste0("Calibration of lambda 46 & dot5: ", set_name),
         subtitle = "(LHS 1000 points; treatment: PE; red: emulator solution)") +
    theme_minimal()
}

# Define population and calibration design
# (shared by all parameter sets, so their calibrations are directly comparable)
pop <- make_population(seed = 5678)
treatment <- list("PE" = c("prrt", "everolimus", NA))
calibration_grid <- make_calibration_design(n_design = 1000, seed = 20260728)

# Parameter sets to calibrate; names are used in the output file names
param_files <- c(
  "NETTER&RADIANT" = "Data/input_params_NETTER&RADIANT.csv",
  "COMPETE" = "Data/input_params_COMPETE.csv"
)

calibrated_params <- list()

for (set_name in names(param_files)) {
  message("=== Parameter set: ", set_name, " (", param_files[[set_name]], ")")
  param_input <- data.table::as.data.table(read.csv(param_files[[set_name]]))
  
  # Calibration results
  calibration_results <- calibrate_lambda_46_and_dot5(
    calibration_grid = calibration_grid,
    param_input = param_input,
    seed = 5678,
    population = pop,
    treatments = treatment,
    n_cores = max(1L, parallel::detectCores() - 1L)
  )
  
  # Relative differences put both targets on the same scale
  calibration_results[
    ,
    `:=`(
      os_median_rel_diff = (os_median - os_median_target) / os_median_target,
      mds_aml_5y_cir_rel_diff = (mds_aml_5y_cir - mds_aml_5y_cir_target) / mds_aml_5y_cir_target
    )
  ][, distance := sqrt(os_median_rel_diff^2 + mds_aml_5y_cir_rel_diff^2)]
  setorder(calibration_results, distance)
  
  write.csv(calibration_results, 
            file.path("Data", paste0("calibrate_lambda_46&dot5_results_", set_name, ".csv"))
            )
  
  emulator_fit <- fit_calibration_emulators(calibration_results, set_name)
  print(emulator_fit$solutions)
  
  write.csv(emulator_fit$solutions,
            file.path("Data", paste0("calibrate_lambda_46&dot5_emulator_", set_name, ".csv")),
            row.names = FALSE)
  
  pdf(file.path("Fig", paste0("Calibrate_lambda_46&dot5_results_", set_name, ".pdf")), width = 5, height = 4)
  print(plot_calibration(calibration_results, emulator_fit$solutions, set_name))
  dev.off()
  
  calibrated_params[[set_name]] <- emulator_fit$calibrated_params
}

# Calibrated est and se to enter in each input_params file
calibrated_params <- rbindlist(calibrated_params, idcol = "parameter_set")
print(calibrated_params)
write.csv(calibrated_params,
          file.path("Data", "calibrate_lambda_46&dot5_input_params.csv"),
          row.names = FALSE)
