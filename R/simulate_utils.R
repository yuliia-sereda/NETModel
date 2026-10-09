#' Sample transition-rate parameters
#'
#' @description
#' Samples person-level transition-rate parameters from log-normal
#' distributions. For median survival times, the sampled time is converted to
#' an exponential rate using `log(2) / time`. For incidence rates, the sampled
#' value is returned directly as the rate.
#'
#' @param N Positive integer giving the number of parameter values to sample.
#' @param est Positive numeric scalar giving the point estimate. This is a
#'   median survival time (years) when `point_stat = "median"` and an incidence rate
#'   (1/year)  when `point_stat = "ir"`.
#' @param se Nonnegative numeric scalar giving the standard error. For median
#'   survival times, it is assumed to be on the original time scale and is
#'   converted approximately to the log scale as `se / est`. For incidence
#'   rates, it is assumed to be on the log-rate scale.
#' @param point_stat Character scalar specifying the type of point estimate.
#'   Supported values are `"median"` and `"ir"`.
#'
#' @return A numeric vector of sampled transition rates of length `N`.
#'
#' @export
sample_lambda <- function(N, est, se, point_stat) {
  if (length(est) != 1L || !is.finite(est) || est <= 0) {
    stop("`est` must be a positive finite number.")
  }
  
  if (length(se) != 1L || !is.finite(se) || se < 0) {
    stop("`se` must be a nonnegative finite number.")
  }
  
  if (!point_stat %in% c("median", "ir")) {
    stop("Unknown `point_stat`: ", point_stat, ". Expected `median` or `ir`.")
  }
  
  if (point_stat == "median") {
    sampled_time <- stats::rlnorm(n = N, meanlog = log(est), sdlog = se / est)
    return(log(2) / sampled_time)
  }
  
  stats::rlnorm(n = N, meanlog = log(est), sdlog = se)
}

#' Sample person-level model parameters
#'
#' @description
#' Samples one set of transition-rate parameters for each simulated person
#' using the estimates and standard errors supplied in `param_input`.
#' Line-specific parameters that are not supplied directly are initialized
#' using equality constraints. Everolimus and PRRT rates are initially assumed
#' to be equal across treatment lines, while the supplied second-line
#' cabozantinib rate is initially assigned to the first and third lines.
#'
#' @param N Positive integer giving the simulated population size.
#' @param param_input A `data.table` containing one row per base parameter and
#'   the columns `param`, `est`, `se`, and `point_stat`.
#'
#' @return A named list of numeric parameter vectors. Each vector has length
#'   `N`, with one sampled value per simulated person.
#'
#' @export
sample_parameters <- function(N, param_input) {
  # Column names used in data.table expressions (declared for R CMD check)
  est <- se <- point_stat <- NULL
  
  stopifnot(all(c("param", "est", "se", "point_stat") %in% colnames(param_input)))
  params <- setNames(
    lapply(seq_len(nrow(param_input)), function(i) {
      sample_lambda(
        N = N,
        est = param_input[i, est],
        se = param_input[i, se],
        point_stat = param_input[i, point_stat]
      )
    }),
    param_input$param
  )
  # Constraints
  params[["lambda_everolimus_34"]] <- params[["lambda_everolimus_23"]] <- params[["lambda_everolimus_12"]]
  params[["lambda_prrt_34"]] <- params[["lambda_prrt_23"]] <- params[["lambda_prrt_12"]]
  params[["lambda_cabozantinib_12"]] <- params[["lambda_cabozantinib_23"]] 
  params[["lambda_cabozantinib_34"]] <- params[["lambda_cabozantinib_23"]] 
  return(params)
}

#' Apply design-specific parameter modifications
#'
#' @description
#' Updates person-level transition-rate parameters for a simulation design
#' point. The function applies treatment-effect discounting to second- and
#' third-line progression rates and multiplies the treatment-related MDS/AML
#' incidence rate by a specified factor.
#'
#' Everolimus and PRRT parameters are supplied as first-line rates. Their
#' second- and third-line rates are obtained by multiplying the first-line rate
#' by the corresponding discounting factor. Cabozantinib is supplied as a
#' second-line rate, so its first- and third-line rates are derived relative to
#' that reference.
#'
#' @param params A named list of person-level model parameter vectors.
#' @param discounting A positive numeric vector of length two. The first value
#'   modifies second-line treatment effects and the second value modifies
#'   third-line treatment effects. The default `c(1, 1)` applies no
#'   discounting.
#' @param aml_mds_multiplier Positive numeric scalar used to multiply
#'   `lambda_dot5`, the treatment-related MDS/AML incidence-rate parameter.
#'   The default value of `1` leaves the rate unchanged.
#'
#' @return The updated named list of model parameter vectors.
#'
#' @export
update_parameters <- function(params, discounting = c(1, 1), aml_mds_multiplier = 1) {
  params[["lambda_everolimus_23"]] <- params[["lambda_everolimus_12"]] * discounting[1]
  params[["lambda_everolimus_34"]] <- params[["lambda_everolimus_12"]] * discounting[2]
  params[["lambda_prrt_23"]] <- params[["lambda_prrt_12"]] * discounting[1]
  params[["lambda_prrt_34"]] <- params[["lambda_prrt_12"]] * discounting[2]
  params[["lambda_cabozantinib_12"]] <- params[["lambda_cabozantinib_23"]] / discounting[1]
  params[["lambda_cabozantinib_34"]] <- (params[["lambda_cabozantinib_23"]] / discounting[1]) * discounting[2]
  params[["lambda_dot5"]] <- params[["lambda_dot5"]] * aml_mds_multiplier
  return(params)
}

#' Simulate NEC outcomes across parameter-design points
#'
#' @description
#' Simulates the supplied population under all specified treatment sequences
#' for each row of the uncertainty design. Person-level base parameters are
#' sampled once and reused across design rows. For each design point, the
#' function applies the corresponding treatment-effect discounting factors and
#' MDS/AML-rate multiplier, initializes a `Population` object, and extracts
#' overall survival, progression-free survival, and 5-year MDS/AML cumulative
#' incidence summaries by treatment sequence.
#'
#' Design rows use common random numbers: every row starts the event
#' simulation from the same random-number state, so differences between rows
#' reflect the design parameters rather than Monte Carlo noise. Results are
#' therefore identical for any value of `n_cores`.
#'
#' @param param_input A `data.table` containing base parameter names, point
#'   estimates, standard errors, and point-statistic types. Required columns
#'   include `param`, `est`, `se`, and `point_stat`.
#' @param design A `data.table` containing one row per simulation design point
#'   and the columns `discount1`, `discount2`, and `aml_mds_multiplier`.
#' @param population A `data.table` describing the simulated population. It
#'   must contain at least `id`, `sex`, `birth_cohort`, `spawn_age`, and
#'   `max_simulation_age`.
#' @param treatments A uniquely named list of treatment sequences accepted by
#'   the `Population` class. The list names are used as treatment-sequence
#'   identifiers in the returned results.
#' @param seed Integer random-number seed used when sampling the person-level
#'   base parameters and the event times.
#' @param show_progress Logical. If `TRUE`, display a text progress bar over
#'   the rows of `design`. Ignored when `n_cores > 1`.
#' @param n_cores Positive integer giving the number of forked worker
#'   processes used to simulate design rows in parallel (via
#'   [parallel::mclapply()], so not available on Windows).
#'
#' @return A `data.table` with one row per design point and treatment sequence.
#'   The output includes mean and median overall survival, mean and median
#'   progression-free survival, associated standard errors, 5-year MDS/AML
#'   cumulative incidence and its standard error, and the design parameters.
#'
#' @export
sim_data <- function(param_input,
                     design,
                     population,
                     treatments,
                     seed,
                     show_progress = TRUE,
                     n_cores = 1L
                     ) {
  stopifnot(all(c("discount1", "discount2", "aml_mds_multiplier") %in% names(design)))
  stopifnot(length(n_cores) == 1L, n_cores >= 1)
  
  N <- nrow(population)
  # Sample parameters
  set.seed(seed)
  param_set <- sample_parameters(N = N, param_input = param_input)
  # Common random numbers: every design row starts from this state.
  rng_state <- get(".Random.seed", envir = globalenv())
  
  simulate_design_row <- function(i) {
    assign(".Random.seed", rng_state, envir = globalenv())
    simulate_design_point(param_set, design, i, population, treatments)
  }
  
  if (n_cores > 1) {
    results <- parallel::mclapply(
      seq_len(nrow(design)),
      simulate_design_row,
      mc.cores = n_cores,
      mc.preschedule = FALSE
    )
    check_parallel_results(results, "Design row")
  } else {
    if (isTRUE(show_progress)) {
      pb <- txtProgressBar(min = 0, max = nrow(design), style = 3)
      on.exit(close(pb), add = TRUE)
    }
    results <- vector("list", nrow(design))
    for (i in seq_len(nrow(design))) {
      results[[i]] <- simulate_design_row(i)
      if (isTRUE(show_progress)) { setTxtProgressBar(pb, i) }
    }
  }
  
  return(rbindlist(results, use.names = TRUE, fill = TRUE))
}

# Simulate one design point and summarise outcomes by treatment sequence.
# `params` holds person-level rate vectors (one value per person); `summarise`
# turns the simulated Population into one row per treatment sequence.
simulate_design_point <- function(params, design, i, population, treatments,
                                  summarise = summarise_population) {
  # Column names used in data.table expressions (declared for R CMD check)
  discount1 <- discount2 <- aml_mds_multiplier <- NULL
  
  params_i <- update_parameters(
    params = params,
    discounting = c(design[i, discount1], design[i, discount2]),
    aml_mds_multiplier = design[i, aml_mds_multiplier]
  )
  
  P <- Population$new(
    population = population,
    treatments = treatments,
    parameters = params_i,
    simulate_at_initialization = TRUE
  )
  
  result <- summarise(P)
  result[, `:=`(
    design_id = i,
    discount1 = design[i, discount1],
    discount2 = design[i, discount2],
    aml_mds_multiplier = design[i, aml_mds_multiplier]
  )]
  result
}

#' Summarise a simulated population: mean and median OS and PFS, and 5-year
#' MDS/AML cumulative incidence, with standard errors
#'
#' @param P A simulated `Population` object.
#'
#' @return A `data.table` with one row per treatment sequence.
#'
#' @export
summarise_population <- function(P) {
  # Column names used in data.table expressions (declared for R CMD check)
  tx <- rmean <- incidence <- se <- `0.95LCL` <- `0.95UCL` <- `se(rmean)` <- NULL
  
  os <- P$getOS()[["surv_time"]][
    , list(
      tx, 
      os_median = median,
      os_median_se = (`0.95UCL` - `0.95LCL`)/3.92,
      os_mean = rmean,
      os_mean_se = `se(rmean)`
    )
  ]
  
  pfs <- P$getPFS()[["surv_time"]][
    , list(
      tx, 
      pfs_median = median,
      pfs_median_se = (`0.95UCL` - `0.95LCL`)/3.92,
      pfs_mean = rmean,
      pfs_mean_se = `se(rmean)`
    )
  ]
  
  mds_aml_5y_cir <- P$getCumIncidenceMDSAML(time_points = 5)[
    , list(
      tx, 
      mds_aml_5y_cir = incidence,
      mds_aml_5y_cir_se = se
    )
  ]
  
  os[pfs, on = "tx"][mds_aml_5y_cir, on = "tx"]
}

#' Summarise a simulated population for the clinician app
#'
#' @description
#' Mean and median OS and PFS (years), OS and PFS survival probabilities at 1,
#' 3 and 5 years, and MDS/AML cumulative incidence at 1, 3 and 5 years (death
#' as a competing event). Probabilities are proportions (0-1).
#'
#' @param P A simulated `Population` object.
#' @param time_points Numeric vector of years for the survival probabilities
#'   and cumulative incidences.
#'
#' @return A `data.table` with one row per treatment sequence.
#'
#' @export
summarise_population_app <- function(P, time_points = c(1, 3, 5)) {
  # Column names used in data.table expressions (declared for R CMD check)
  tx <- rmean <- NULL
  
  survival_summary <- function(km, prefix) {
    times <- km[["surv_time"]][, list(tx, mean = rmean, median = median)]
    setnames(times, c("mean", "median"), paste0(prefix, c("_mean", "_median")))
    probs <- dcast(km[["surv_prop"]], tx ~ time, value.var = "survival")
    setnames(probs, setdiff(names(probs), "tx"), paste0(prefix, "_", setdiff(names(probs), "tx"), "y"))
    times[probs, on = "tx"]
  }
  os <- survival_summary(P$getOS(time_points = time_points), "os")
  pfs <- survival_summary(P$getPFS(time_points = time_points), "pfs")
  cir <- dcast(P$getCumIncidenceMDSAML(time_points = time_points), tx ~ time, value.var = "incidence")
  setnames(cir, setdiff(names(cir), "tx"), paste0("mds_aml_", setdiff(names(cir), "tx"), "y_cir"))
  os[pfs, on = "tx"][cir, on = "tx"]
}

# A parallel worker that errors returns a try-error; one that is killed (e.g.
# out of memory) returns NULL. Treat both as failures rather than dropping rows.
# `ids` labels the results in error messages (e.g. parameter-set IDs).
check_parallel_results <- function(results, what, ids = seq_along(results)) {
  failed <- vapply(results, function(x) is.null(x) || inherits(x, "try-error"), logical(1))
  if (any(failed)) {
    stop(what, "(s) ", paste(ids[failed], collapse = ", "), " failed: ",
         results[[which(failed)[1]]])
  }
  invisible(TRUE)
}

#' Simulate outcomes under parameter uncertainty (probabilistic sensitivity analysis)
#'
#' @description
#' Draws `n_sets` sets of model parameters from their uncertainty
#' distributions (one value per parameter per set) and, for each set,
#' simulates the whole population with those parameter values at every
#' design point. Outcomes are summarised by set, treatment sequence, and
#' design point, so their spread across sets reflects parameter uncertainty
#' rather than Monte Carlo error.
#'
#' Each set has its own random-number stream; within a set, all design points
#' start from the same state (common random numbers). Results are therefore
#' identical for any value of `n_cores`.
#'
#' @inheritParams sim_data
#' @param n_sets Positive integer giving the number of parameter sets.
#' @param seed Integer random-number seed for the parameter draws and the
#'   per-set random-number streams.
#' @param n_cores Positive integer giving the number of forked worker
#'   processes used to simulate parameter sets in parallel (via
#'   [parallel::mclapply()], so not available on Windows).
#' @param set_ids Integer vector of the parameter sets (from `1:n_sets`) to
#'   simulate, e.g. one chunk per array job. All `n_sets` sets are always drawn
#'   from `seed`, so a set gives the same result whichever chunk it is in.
#' @param summarise Function turning a simulated `Population` into one row per
#'   treatment sequence, e.g. [summarise_population()] (default) or
#'   [summarise_population_app()].
#'
#' @return A `data.table` with one row per parameter set, design point, and
#'   treatment sequence. It contains `set_id`, the outcome summaries returned
#'   by `summarise`, the design parameters, and the sampled base parameter
#'   values of the set (before design-specific adjustment).
#'
#' @export
sim_parameter_sets <- function(param_input,
                               design,
                               population,
                               treatments,
                               n_sets,
                               seed,
                               n_cores = 1L,
                               set_ids = seq_len(n_sets),
                               summarise = summarise_population) {
  # Column names used in data.table expressions (declared for R CMD check)
  set_id <- NULL
  
  stopifnot(all(c("discount1", "discount2", "aml_mds_multiplier") %in% names(design)))
  stopifnot(length(n_sets) == 1L, n_sets >= 1)
  stopifnot(length(n_cores) == 1L, n_cores >= 1)
  stopifnot(length(set_ids) >= 1L, all(set_ids %in% seq_len(n_sets)), !anyDuplicated(set_ids))
  
  N <- nrow(population)
  set.seed(seed)
  # One value per parameter and set
  parameter_draws <- sample_parameters(N = n_sets, param_input = param_input)
  set_seeds <- sample.int(.Machine$integer.max, n_sets)
  
  simulate_set <- function(k) {
    # Everyone in the population shares the set's parameter values
    params_k <- lapply(parameter_draws, function(v) rep(v[k], N))
    set.seed(set_seeds[k])
    rng_state <- get(".Random.seed", envir = globalenv())
    result <- rbindlist(lapply(seq_len(nrow(design)), function(i) {
      assign(".Random.seed", rng_state, envir = globalenv())
      simulate_design_point(params_k, design, i, population, treatments, summarise)
    }))
    result[, set_id := k]
    result
  }
  
  if (n_cores > 1) {
    results <- parallel::mclapply(set_ids, simulate_set,
                                  mc.cores = n_cores, mc.preschedule = FALSE)
    check_parallel_results(results, "Parameter set", ids = set_ids)
  } else {
    results <- lapply(set_ids, simulate_set)
  }
  
  base_params <- param_input$param
  draws <- as.data.table(parameter_draws[base_params])
  draws[, set_id := seq_len(n_sets)]
  draws <- draws[set_id %in% set_ids]
  
  out <- rbindlist(results, use.names = TRUE)
  out <- draws[out, on = "set_id"]
  setcolorder(out, c("set_id", setdiff(names(out), c("set_id", base_params)), base_params))
  out[]
}

#' Treatment sequences compared in the analysis
#'
#' @description
#' Returns the eight treatment sequences used throughout the analysis, as a
#' uniquely named list accepted by the `Population` class. Names are the
#' sequence codes (P = PRRT, E = everolimus, C = cabozantinib); two-drug
#' sequences have `NA` as the third treatment.
#'
#' @return A named list of character vectors of length three.
#'
#' @export
treatment_sequences <- function() {
  list(
    # Cabozantinib is not available
    "PE" = c("prrt", "everolimus", NA),
    "EP" = c("everolimus", "prrt", NA),
    # Start with PRRT, Cabozantinib is available
    "PEC" = c("prrt", "everolimus", "cabozantinib"),
    "PCE" = c("prrt", "cabozantinib", "everolimus"),
    # Start with Everolimus, Cabozantinib is available
    "EPC" = c("everolimus", "prrt", "cabozantinib"),
    "ECP" = c("everolimus", "cabozantinib", "prrt"),
    # Start with Cabozantinib
    "CPE" = c("cabozantinib", "prrt", "everolimus"),
    "CEP" = c("cabozantinib", "everolimus", "prrt")
  )
}

#' Create the simulated NEC population
#'
#' @description
#' Creates the population of `N` people used in the simulations: sex is
#' sampled 1:1, all people belong to the 1984 birth cohort, and model-entry
#' age is drawn from a normal distribution with mean 63 and SD 9, truncated
#' to 18-84 years (NETTER-1).
#'
#' @param N Positive integer giving the population size (default 100,000).
#' @param seed Integer random-number seed.
#' @param max_simulation_age Numeric maximum simulated age, in years.
#'
#' @return A `data.table` with the columns `id`, `sex`, `birth_cohort`,
#'   `spawn_age`, and `max_simulation_age`.
#'
#' @export
make_population <- function(N = 100000, seed, max_simulation_age = 100) {
  set.seed(seed)
  data.table::data.table(
    id = seq_len(N),
    sex = sample(x = c("male", "female"), size = N, prob = c(0.5, 0.5), replace = TRUE),
    birth_cohort = rep(1984, N),
    spawn_age = truncnorm::rtruncnorm(n = N, mean = 63, sd = 9, a = 18, b = 84),
    max_simulation_age = max_simulation_age
  )
}
