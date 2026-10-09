# Input parameters on the COMPETE scale, with strong attenuation so that the
# thinning majorizers differ a lot between treatments in later lines.
regression_param_input <- function(lambda_dot5 = 0.0039) {
  data.table::data.table(
    param = c("lambda_everolimus_12", "lambda_prrt_12", "lambda_cabozantinib_23",
              "lambda_dot5", "lambda_57", "lambda_46"),
    point_stat = c("median", "median", "median", "ir", "median", "median"),
    est = c(1, 1.991666667, 0.7, lambda_dot5, 1, 0.62),
    se = c(0.316751701, 0.316751701, 0.108418367, 0.278413, 0.331632653, 0.0425)
  )
}

simulate_regression_population <- function(N = 20000, lambda_dot5 = 0.0039,
                                           discounting = c(1.5, 2.25), aml_mds_multiplier = 1) {
  pop <- make_population(N = N, seed = 1)
  set.seed(2)
  pars <- update_parameters(
    sample_parameters(N, regression_param_input(lambda_dot5)),
    discounting = discounting,
    aml_mds_multiplier = aml_mds_multiplier
  )
  set.seed(3)
  Population$new(population = pop, parameters = pars)
}

# Rate of the treatment given in each row of `tx` for a line ("12", "23", "34")
line_rate <- function(d, tx, line) {
  data.table::fcase(
    tx == "prrt", d[[paste0("params_lambda_prrt_", line)]],
    tx == "everolimus", d[[paste0("params_lambda_everolimus_", line)]],
    tx == "cabozantinib", d[[paste0("params_lambda_cabozantinib_", line)]]
  )
}

# Simulated share of units without an event vs. its analytic expectation,
# E[exp(-rate * time at risk)]; returns the difference in standard errors.
share_z <- function(no_event, rate, time_at_risk) {
  p_sim <- mean(no_event)
  p_exp <- mean(exp(-rate * time_at_risk))
  (p_sim - p_exp) / sqrt(p_exp * (1 - p_exp) / length(no_event))
}

# With a negligible MDS/AML risk, every progression is observed, so the share
# of people who do not progress before other-cause death has a closed form.
P_no_mds <- simulate_regression_population(lambda_dot5 = 1e-10)
d_no_mds <- P_no_mds$getPopulation()

test_that("line 1 progression matches its analytic expectation for every first treatment", {
  d <- unique(d_no_mds, by = c("id", "treatment1"))
  d[, rate := line_rate(d, treatment1, "12")]
  for (tx in unique(d$treatment1)) {
    x <- d[treatment1 == tx]
    z <- share_z(is.na(x$age_progression_to_stable_2), x$rate,
                 x$age_dead_from_other_causes - x$spawn_age)
    expect_lt(abs(z), 4, label = paste("line 1,", tx))
  }
})

test_that("line 2 progression matches its analytic expectation for every second treatment", {
  d <- unique(d_no_mds[!is.na(age_progression_to_stable_2)], by = c("id", "treatment1", "treatment2"))
  d[, rate := line_rate(d, treatment2, "23")]
  # Two-drug sequences move to salvage, three-drug sequences to stable 3.
  d[, no_event := data.table::fifelse(is.na(treatment3),
                                      is.na(age_progression_to_salvage_therapy),
                                      is.na(age_progression_to_stable_3))]
  for (tx in unique(d$treatment2)) {
    x <- d[treatment2 == tx]
    z <- share_z(x$no_event, x$rate, x$age_dead_from_other_causes - x$age_progression_to_stable_2)
    expect_lt(abs(z), 4, label = paste("line 2,", tx))
  }
})

test_that("line 3 progression to salvage matches its analytic expectation (PRRT in line 3 included)", {
  # Before the fix, the line-3 majorizer was the population maximum over all
  # treatments; with atmostB = 20 candidates, PRRT-third-line events were
  # dropped and far too few people reached salvage (12.5% vs 5.5% expected).
  d <- d_no_mds[!is.na(treatment3) & !is.na(age_progression_to_stable_3)]
  d[, rate := line_rate(d, treatment3, "34")]
  for (tx in unique(d$treatment3)) {
    x <- d[treatment3 == tx]
    z <- share_z(is.na(x$age_progression_to_salvage_therapy), x$rate,
                 x$age_dead_from_other_causes - x$age_progression_to_stable_3)
    expect_lt(abs(z), 4, label = paste("line 3,", tx))
  }
})

test_that("death after salvage matches its analytic expectation", {
  d <- d_no_mds[!is.na(age_progression_to_salvage_therapy)]
  z <- share_z(is.na(d$age_dead_from_disease_progression), d$params_lambda_46,
               d$age_dead_from_other_causes - d$age_progression_to_salvage_therapy)
  expect_lt(abs(z), 4)
})

# Realistic MDS/AML risk (unfavourable multiplier) for the remaining tests
P <- simulate_regression_population(aml_mds_multiplier = 3)
d <- P$getPopulation()

test_that("MDS/AML onset matches the dose-capped intensity (doses stop at PRRT stop)", {
  # Only PRRT-first sequences: onset is then not preceded by other events.
  x <- unique(d[treatment1 == "prrt"], by = c("id", "prrt_start_age", "prrt_stop_age"))
  # Cumulative hazard from PRRT start to other-cause death: dose k is given at
  # start + (k - 1) / 6 if PRRT has not stopped, and adds lambda from then on.
  x[, cum_hazard := {
    end <- age_dead_from_other_causes
    stop_age <- pmin(prrt_stop_age, end)
    h <- 0
    for (k in 1:4) {
      given <- prrt_start_age + (k - 1) / 6
      h <- h + ifelse(given < stop_age, pmax(end - given, 0), 0)
    }
    params_lambda_dot5 * h
  }]
  z <- share_z(is.na(x$latent_age_with_mds_aml), rep(1, nrow(x)), x$cum_hazard)
  expect_lt(abs(z), 4)
})

test_that("intensity_dot5 counts only the PRRT doses received before PRRT stops", {
  t <- c(59, 60, 60.1, 60.2, 60.3, 61, 70)
  expect_equal(intensity_dot5(t, 1, prrt_start_age = 60), c(0, 0, 1, 2, 2, 4, 4))
  expect_equal(intensity_dot5(t, 1, prrt_start_age = 60, prrt_stop_age = 60.25), c(0, 0, 1, 2, 2, 2, 2))
  expect_equal(intensity_dot5(t, 1, prrt_start_age = Inf), rep(0, length(t)))
})

test_that("no treatment runs, starts, or progresses after MDS/AML", {
  m <- d[!is.na(age_with_mds_aml)]
  expect_gt(nrow(m), 100)
  for (tx in c("prrt", "everolimus", "cabozantinib")) {
    start <- m[[paste0(tx, "_start_age")]]
    stop <- m[[paste0(tx, "_stop_age")]]
    expect_false(any(is.finite(start) & stop > m$age_with_mds_aml), label = paste(tx, "running"))
    expect_false(any(is.finite(start) & start >= m$age_with_mds_aml), label = paste(tx, "started"))
  }
  for (v in c("age_progression_to_stable_2", "age_progression_to_stable_3",
              "age_progression_to_salvage_therapy", "age_dead_from_disease_progression")) {
    expect_false(any(!is.na(m[[v]]) & m[[v]] >= m$age_with_mds_aml), label = v)
  }
})

test_that("5-year MDS/AML incidence treats death as a competing event", {
  # With no censoring, the competing-risk cumulative incidence is the plain
  # proportion with MDS/AML by 5 years. Treating death as censoring (the old
  # code) gives 1 - Kaplan-Meier, which is larger.
  cif <- P$getCumIncidenceMDSAML(time_points = 5)
  expected <- d[, .(expected = mean(!is.na(age_with_mds_aml) & age_with_mds_aml - spawn_age <= 5)),
                by = .(tx = tx_seq)]
  x <- merge(cif, expected, by = "tx")
  expect_equal(nrow(x), length(treatment_sequences()))
  expect_equal(x$incidence, x$expected, tolerance = 1e-8)
})

test_that("cumulative incidence is not rounded and has a standard error", {
  # The old code rounded the incidence to 3 decimals before it reached the GAM.
  cif <- P$getCumIncidenceMDSAML(time_points = 5)
  expect_true(any(abs(cif$incidence * 1000 - round(cif$incidence * 1000)) > 1e-9))
  expect_true(all(cif$se > 0))
})

test_that("PFS ends at the earliest of salvage, MDS/AML, and death; OS at death", {
  pfs <- P$getPFS()$surv_time[, .(tx, rmean)]
  expected_pfs <- d[, .(expected = mean(pmin(age_progression_to_salvage_therapy, age_with_mds_aml,
                                             age_dead_from_all_causes, na.rm = TRUE) - spawn_age)),
                    by = .(tx = tx_seq)]
  x <- merge(pfs, expected_pfs, by = "tx")
  expect_equal(x$rmean, x$expected, tolerance = 1e-6)

  os <- P$getOS()$surv_time[, .(tx, rmean)]
  expected_os <- d[, .(expected = mean(age_dead_from_all_causes - spawn_age)), by = .(tx = tx_seq)]
  x <- merge(os, expected_os, by = "tx")
  expect_equal(x$rmean, x$expected, tolerance = 1e-6)
})

test_that("survival summaries keep non-integer time points", {
  # The old code converted time points with as.integer(), so 1.5 became 1.
  expect_true(1.5 %in% P$getPFS(time_points = c(1, 1.5))$surv_prop$time)
})

test_that("reset() allows the population to be simulated again", {
  N <- 500
  pop <- make_population(N = N, seed = 4)
  set.seed(5)
  pars <- sample_parameters(N, regression_param_input())
  P_small <- Population$new(population = pop, parameters = pars)
  expect_no_error(P_small$reset()$update())
  d_small <- P_small$getPopulation()
  expect_equal(nrow(d_small), N * length(treatment_sequences()))
  expect_equal(data.table::uniqueN(d_small, by = c("id", "tx_seq")), nrow(d_small))
})

test_that("make_population() matches the paper's population", {
  expect_equal(formals(make_population)$N, 100000)
  pop <- make_population(N = 20000, seed = 1)
  expect_true(all(pop$spawn_age >= 18 & pop$spawn_age <= 84))
  expect_equal(mean(pop$spawn_age), 63, tolerance = 0.01)
  expect_equal(mean(pop$sex == "female"), 0.5, tolerance = 0.03)
  expect_true(all(pop$max_simulation_age == 100))
  expect_named(treatment_sequences(), c("PE", "EP", "PEC", "PCE", "EPC", "ECP", "CPE", "CEP"))
})

test_that("sample_lambda() uses se / est for medians and se for log rates", {
  set.seed(6)
  x <- sample_lambda(2e5, est = 2, se = 0.4, point_stat = "median")
  expect_equal(mean(log(x)), log(log(2)) - log(2), tolerance = 0.01)
  expect_equal(sd(log(x)), 0.4 / 2, tolerance = 0.01)
  y <- sample_lambda(2e5, est = 0.004, se = 0.28, point_stat = "ir")
  expect_equal(mean(log(y)), log(0.004), tolerance = 0.01)
  expect_equal(sd(log(y)), 0.28, tolerance = 0.01)
})

test_that("update_parameters() applies cumulative position factors", {
  params <- list(
    lambda_everolimus_12 = 1, lambda_prrt_12 = 0.5, lambda_cabozantinib_23 = 2, lambda_dot5 = 0.01
  )
  u <- update_parameters(params, discounting = c(1.5, 2.25), aml_mds_multiplier = 3)
  expect_equal(u$lambda_everolimus_23, 1.5)
  expect_equal(u$lambda_everolimus_34, 2.25)
  expect_equal(u$lambda_prrt_34, 0.5 * 2.25)
  # Cabozantinib is observed in position 2
  expect_equal(u$lambda_cabozantinib_12, 2 / 1.5)
  expect_equal(u$lambda_cabozantinib_34, 2 * 2.25 / 1.5)
  expect_equal(u$lambda_dot5, 0.03)
})

test_that("sim_data() uses outcome-first names and common random numbers", {
  N <- 2000
  pop <- make_population(N = N, seed = 7)
  design <- data.table::data.table(discount1 = c(1, 1.5, 1), discount2 = c(1, 2.25, 1),
                                   aml_mds_multiplier = c(1, 3, 1))
  res <- sim_data(regression_param_input(), design, pop, treatment_sequences(), seed = 8,
                  show_progress = FALSE)
  expect_true(all(c("os_mean", "os_median", "pfs_mean", "pfs_median", "mds_aml_5y_cir",
                    "os_mean_se", "os_median_se", "pfs_mean_se", "pfs_median_se",
                    "mds_aml_5y_cir_se") %in% names(res)))
  expect_false(any(grepl("^(mean|median)_|mdsaml", names(res))))
  expect_equal(nrow(res), nrow(design) * length(treatment_sequences()))
  # Identical design rows give identical results
  cols <- setdiff(names(res), "design_id")
  expect_equal(res[design_id == 1, ..cols], res[design_id == 3, ..cols])

  skip_on_os("windows")
  skip_on_cran()
  res2 <- sim_data(regression_param_input(), design, pop, treatment_sequences(), seed = 8,
                   show_progress = FALSE, n_cores = 2)
  expect_equal(res2, res)
})

test_that("sim_parameter_sets() propagates parameter uncertainty, not Monte Carlo error", {
  # The bounds simulation used to repeat sim_data() with different seeds. Rates
  # were drawn per person, so each run averaged over parameter uncertainty and
  # the spread across runs was only Monte Carlo error.
  N <- 2000
  n_sets <- 6
  pop <- make_population(N = N, seed = 9)
  design <- data.table::data.table(discount1 = c(1, 1.5, 1), discount2 = c(1, 2.25, 1),
                                   aml_mds_multiplier = c(1, 3, 1))
  res <- sim_parameter_sets(regression_param_input(), design, pop, treatment_sequences(),
                            n_sets = n_sets, seed = 10)
  expect_equal(nrow(res), n_sets * nrow(design) * length(treatment_sequences()))
  base_params <- regression_param_input()$param
  expect_true(all(c("set_id", "os_mean", "pfs_mean", "mds_aml_5y_cir", base_params) %in% names(res)))
  # One value of each parameter per set, differing between sets
  for (p in base_params) {
    expect_true(all(res[, data.table::uniqueN(get(p)), by = set_id]$V1 == 1), label = p)
  }
  expect_equal(data.table::uniqueN(res$lambda_prrt_12), n_sets)
  # Within a set, identical design points share random numbers
  cols <- setdiff(names(res), "design_id")
  expect_equal(res[design_id == 1, ..cols], res[design_id == 3, ..cols])
  # Parameter uncertainty dominates Monte Carlo error across sets
  mc_runs <- data.table::rbindlist(lapply(1:n_sets, function(i) {
    sim_data(regression_param_input(), design[1], pop, treatment_sequences(), seed = i,
             show_progress = FALSE)
  }))
  expect_gt(res[tx == "PE" & design_id == 1, sd(os_mean)], 3 * mc_runs[tx == "PE", sd(os_mean)])

  skip_on_os("windows")
  skip_on_cran()
  res2 <- sim_parameter_sets(regression_param_input(), design, pop, treatment_sequences(),
                             n_sets = n_sets, seed = 10, n_cores = 2)
  expect_equal(res2, res)
})

test_that("summarise_population_app() matches the simulated event times", {
  s <- summarise_population_app(P)
  expect_named(s, c("tx", app_outcome_info()$outcome), ignore.order = TRUE)
  expect_equal(nrow(s), length(treatment_sequences()))
  expected <- d[, {
    os_time <- age_dead_from_all_causes - spawn_age
    pfs_time <- pmin(age_progression_to_salvage_therapy, age_with_mds_aml,
                     age_dead_from_all_causes, na.rm = TRUE) - spawn_age
    mds_time <- age_with_mds_aml - spawn_age
    .(os_mean = mean(os_time), os_3y = mean(os_time > 3), pfs_5y = mean(pfs_time > 5),
      mds_aml_1y_cir = mean(!is.na(mds_time) & mds_time <= 1))
  }, by = .(tx = tx_seq)]
  x <- merge(s, expected, by = "tx", suffixes = c("", ".expected"))
  for (v in c("os_mean", "os_3y", "pfs_5y", "mds_aml_1y_cir")) {
    expect_equal(x[[v]], x[[paste0(v, ".expected")]], tolerance = 1e-6, label = v)
  }
})

test_that("sim_parameter_sets() gives the same results when run in chunks of sets", {
  N <- 1000
  pop <- make_population(N = N, seed = 11)
  design <- data.table::data.table(discount1 = c(1, 1.2), discount2 = c(1, 1.56),
                                   aml_mds_multiplier = c(1, 2))
  run <- function(set_ids) {
    sim_parameter_sets(regression_param_input(), design, pop, treatment_sequences(), n_sets = 5,
                       seed = 12, set_ids = set_ids, summarise = summarise_population_app)
  }
  full <- run(1:5)
  expect_equal(rbind(run(1:2), run(3:5)), full)
  expect_error(run(c(1, 6)))
})

test_that("app summaries report units, 95% intervals and the right direction of benefit", {
  N <- 1000
  pop <- make_population(N = N, seed = 13)
  design <- data.table::data.table(design_id = 1:2, discount1 = c(1, 1.5), discount2 = c(1, 2.25),
                                   aml_mds_multiplier = c(1, 3))
  sets <- sim_parameter_sets(regression_param_input(), design[, .(discount1, discount2, aml_mds_multiplier)],
                             pop, treatment_sequences(), n_sets = 20, seed = 14,
                             summarise = summarise_population_app)
  info <- app_outcome_info()
  
  o <- summarise_app_outcomes(sets, design = design)
  expect_equal(nrow(o), nrow(design) * length(treatment_sequences()) * nrow(info))
  expect_true(all(o$lower95 <= o$mean & o$mean <= o$upper95))
  expect_true(all(o[unit == "%", mean >= 0 & mean <= 100]))
  pec_os <- sets[tx == "PEC" & design_id == 1, os_mean]
  expect_equal(o[tx == "PEC" & design_id == 1 & outcome == "os_mean", mean], mean(pec_os))
  expect_equal(o[tx == "PEC" & design_id == 1 & outcome == "os_mean", upper95],
               unname(quantile(pec_os, 0.975)))
  
  df <- summarise_app_differences(sets, design = design)
  expect_equal(nrow(df), nrow(design) * choose(length(treatment_sequences()), 2) * nrow(info))
  # PEC vs CEP, 5-year MDS/AML CIR: difference in percentage points, and PEC is
  # better when its incidence is lower
  x <- data.table::dcast(sets[design_id == 2 & tx %in% c("PEC", "CEP")], set_id ~ tx,
                         value.var = "mds_aml_5y_cir")
  r <- df[design_id == 2 & tx1 == "PEC" & tx2 == "CEP" & outcome == "mds_aml_5y_cir"]
  expect_equal(r$mean_diff, mean((x$PEC - x$CEP) * 100))
  expect_equal(r$p_tx1_better, mean(x$PEC < x$CEP))
  expect_equal(r$diff_unit, "percentage points")
  expect_equal(df[outcome == "os_mean", unique(diff_unit)], "months")
})
