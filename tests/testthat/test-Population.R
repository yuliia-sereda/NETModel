# Input parameters
set.seed(12345)

N <- 1000L

input_pop <- data.table::data.table(
  id = seq_len(N),
  sex = sample(
    c("male", "female"),
    size = N,
    prob = c(0.5, 0.5),
    replace = TRUE
  ),
  birth_cohort = 1984,
  spawn_age = 60,
  max_simulation_age = 110
)

original_pop <- data.table::copy(input_pop)

parameters <- list(
  lambda_everolimus_12 = 1 / (11 / 12),
  lambda_everolimus_23 = 1 / (11 / 12),
  lambda_everolimus_34 = 1 / (11 / 12),
  lambda_prrt_12 = 1 / (28.4 / 12),
  lambda_prrt_23 = 1 / (28.4 / 12),
  lambda_prrt_34 = 1 / (28.4 / 12),
  lambda_cabozantinib_12 = 1 / (8.4 / 12),
  lambda_cabozantinib_23 = 1 / (8.4 / 12),
  lambda_cabozantinib_34 = 1 / (8.4 / 12),
  lambda_46 = 1 / (6 / 12),
  lambda_57 = 1 / (6 / 12),
  lambda_dot5 = (0.002564665 * 6) / 4
)

single_treatment_sequence <- list(
  PE = c("prrt", "everolimus", NA)
)

multiple_treatment_sequences <- list(
  PE = c("prrt", "everolimus", NA),
  PEC = c("prrt", "everolimus", "cabozantinib"),
  PCE = c("prrt", "cabozantinib", "everolimus")
)

test_that("Population initializes, updates, and protects its input", {
  
  expect_no_error({
    P <- Population$new(
      population = input_pop,
      treatments = multiple_treatment_sequences,
      parameters = parameters,
      simulate_at_initialization = FALSE
    )
    
    P$update()
  })
  
  expect_error(P$update())
  
  # The input object must remain unchanged.
  expect_identical(input_pop, original_pop)
  
  # Automatic initialization should work.
  set.seed(12345)
  
  expect_no_error({
    P1 <- Population$new(
      population = input_pop,
      treatments = single_treatment_sequence,
      parameters = parameters,
      simulate_at_initialization = TRUE
    )
  })
  
  expect_error(P1$update())
})

test_that("Population is expanded to one row per person and treatment sequence", {
  P <- Population$new(
    population = input_pop,
    treatments = multiple_treatment_sequences,
    parameters = parameters,
    simulate_at_initialization = TRUE
  )
  simulated <- P$getPopulation()
  expect_equal(nrow(simulated), N * length(multiple_treatment_sequences))
  expect_equal(data.table::uniqueN(simulated$id), N)
  expect_setequal(unique(simulated$tx_seq), names(multiple_treatment_sequences))
  
  expect_true(
    simulated[
      ,
      all(.N == length(multiple_treatment_sequences)),
      by = id
    ][, all(V1)]
  )
  
  expect_false(
    anyDuplicated(simulated[, .(id, tx_seq)]) > 0L
  )
})

test_that("Baseline quantities are shared across counterfactual strategies", {
  P <- Population$new(
    population = input_pop,
    treatments = multiple_treatment_sequences,
    parameters = parameters,
    simulate_at_initialization = TRUE
  )
  
  simulated <- P$getPopulation()
  
  shared_vars <- c(
    "sex",
    "birth_cohort",
    "spawn_age",
    "max_simulation_age",
    "age_dead_from_other_causes",
    "params_lambda_everolimus_12",
    "params_lambda_everolimus_23",
    "params_lambda_everolimus_34",
    "params_lambda_prrt_12",
    "params_lambda_prrt_23",
    "params_lambda_prrt_34",
    "params_lambda_cabozantinib_12",
    "params_lambda_cabozantinib_23",
    "params_lambda_cabozantinib_34",
    "params_lambda_46",
    "params_lambda_57",
    "params_lambda_dot5"
  )
  
  checks <- simulated[
    ,
    lapply(.SD, data.table::uniqueN),
    by = id,
    .SDcols = shared_vars
  ]
  
  expect_true(
    checks[
      ,
      all(unlist(.SD) == 1L),
      .SDcols = shared_vars
    ]
  )
})

test_that("Shared first-line progression is identical before later competing events", {
  P <- Population$new(
    population = input_pop,
    treatments = multiple_treatment_sequences,
    parameters = parameters,
    simulate_at_initialization = TRUE
  )
  
  simulated <- P$getPopulation()
  
  checks <- simulated[
    ,
    .(
      n_progression_12 = uniqueN(
        age_progression_to_stable_2,
        na.rm = FALSE
      )
    ),
    by = id
  ]
  
  expect_true(checks[, all(n_progression_12 == 1L)])
})

test_that("Two-treatment sequences skip stable state 3", {
  P <- Population$new(
    population = input_pop,
    treatments = multiple_treatment_sequences,
    parameters = parameters,
    simulate_at_initialization = TRUE
  )
  
  simulated <- P$getPopulation()
  
  pe <- simulated[tx_seq == "PE"]
  
  expect_true(pe[, all(is.na(age_progression_to_stable_3))])
  
  # PE should reach salvage from the event after treatment 2.
  expect_true(
    pe[
      !is.na(age_progression_to_salvage_therapy),
      all(
        everolimus_stop_age ==
          age_progression_to_salvage_therapy
      )
    ]
  )
})

test_that("Three-treatment sequences enter stable state 3", {
  P <- Population$new(
    population = input_pop,
    treatments = multiple_treatment_sequences,
    parameters = parameters,
    simulate_at_initialization = TRUE
  )
  
  simulated <- P$getPopulation()
  
  three_drug <- simulated[!is.na(treatment3)]
  
  observed <- three_drug[
    !is.na(age_progression_to_stable_3)
  ]
  
  expect_gt(nrow(observed), 0L)
  
  expect_true(
    observed[
      ,
      all(
        age_progression_to_stable_3 >=
          age_progression_to_stable_2
      )
    ]
  )
  
  expect_true(
    observed[
      !is.na(age_progression_to_salvage_therapy),
      all(
        age_progression_to_salvage_therapy >=
          age_progression_to_stable_3
      )
    ]
  )
})

test_that("MDS/AML latent age is shared for identical PRRT histories", {
  P <- Population$new(
    population = input_pop,
    treatments = multiple_treatment_sequences,
    parameters = parameters,
    simulate_at_initialization = TRUE
  )
  
  simulated <- P$getPopulation()
  
  checks <- simulated[
    is.finite(prrt_start_age),
    .(
      n_latent_mdsaml = uniqueN(
        latent_age_with_mds_aml,
        na.rm = FALSE
      )
    ),
    by = .(
      id,
      prrt_start_age
    )
  ]
  
  expect_true(checks[, all(n_latent_mdsaml == 1L)])
})

test_that("Simulated event ages are chronologically valid", {
  P <- Population$new(
    population = input_pop,
    treatments = multiple_treatment_sequences,
    parameters = parameters,
    simulate_at_initialization = TRUE
  )
  
  simulated <- P$getPopulation()
  
  expect_true(simulated[, all(age_dead_from_all_causes >= spawn_age)])
  expect_true(simulated[, all(age_dead_from_all_causes <= max_simulation_age)])
  expect_true(simulated[, all(age_dead_from_all_causes <= age_dead_from_other_causes)])
  
  expect_true(
    simulated[
      !is.na(age_progression_to_stable_2),
      all(age_progression_to_stable_2 >= spawn_age)
    ]
  )
  
  expect_true(
    simulated[
      !is.na(age_progression_to_stable_3),
      all(
        age_progression_to_stable_3 >=
          age_progression_to_stable_2
      )
    ]
  )
  
  expect_true(
    simulated[
      !is.na(age_progression_to_salvage_therapy),
      all(
        age_progression_to_salvage_therapy >=
          age_progression_to_stable_2
      )
    ]
  )
  
  expect_true(
    simulated[
      !is.na(age_dead_from_disease_progression),
      all(
        age_dead_from_disease_progression >=
          age_progression_to_salvage_therapy
      )
    ]
  )
  
  expect_true(
    simulated[
      !is.na(age_with_mds_aml),
      all(age_with_mds_aml > prrt_start_age)
    ]
  )
  
  expect_true(
    simulated[
      !is.na(age_dead_from_mds_aml),
      all(age_dead_from_mds_aml >= age_with_mds_aml)
    ]
  )
})

test_that("Treatment start and stop ages are consistent", {
  P <- Population$new(
    population = input_pop,
    treatments = multiple_treatment_sequences,
    parameters = parameters,
    simulate_at_initialization = TRUE
  )
  
  simulated <- P$getPopulation()
  
  treatments <- c(
    "prrt",
    "everolimus",
    "cabozantinib"
  )
  
  for (tx in treatments) {
    start_var <- paste0(tx, "_start_age")
    stop_var <- paste0(tx, "_stop_age")
    
    expect_true(
      simulated[
        is.finite(get(start_var)),
        all(is.finite(get(stop_var)))
      ],
      info = paste(tx, "has finite start but non-finite stop")
    )
    
    expect_true(
      simulated[
        is.finite(get(start_var)) &
          is.finite(get(stop_var)),
        all(get(stop_var) >= get(start_var))
      ],
      info = paste(tx, "stop age occurs before start age")
    )
    
    expect_true(
      simulated[
        is.infinite(get(start_var)),
        all(is.infinite(get(stop_var)))
      ],
      info = paste(tx, "never started but has finite stop age")
    )
  }
})

test_that("Treatment-age columns agree with assigned treatment sequence", {
  P <- Population$new(
    population = input_pop,
    treatments = multiple_treatment_sequences,
    parameters = parameters,
    simulate_at_initialization = TRUE
  )
  
  simulated <- P$getPopulation()
  
  expect_true(
    simulated[
      treatment1 == "prrt",
      all(is.finite(prrt_start_age))
    ]
  )
  
  expect_true(
    simulated[
      treatment1 == "everolimus",
      all(is.finite(everolimus_start_age))
    ]
  )
  
  expect_true(
    simulated[
      treatment1 == "cabozantinib",
      all(is.finite(cabozantinib_start_age))
    ]
  )
  
  # Treatments absent from the sequence must never start.
  expect_true(
    simulated[
      !(
        treatment1 == "prrt" |
          treatment2 == "prrt" |
          treatment3 == "prrt"
      ),
      all(is.infinite(prrt_start_age))
    ]
  )
  
  expect_true(
    simulated[
      !(
        treatment1 == "everolimus" |
          treatment2 == "everolimus" |
          treatment3 == "everolimus"
      ),
      all(is.infinite(everolimus_start_age))
    ]
  )
  
  expect_true(
    simulated[
      !(
        treatment1 == "cabozantinib" |
          treatment2 == "cabozantinib" |
          treatment3 == "cabozantinib"
      ),
      all(is.infinite(cabozantinib_start_age))
    ]
  )
})

test_that("Treatment sequences are validated", {
  expect_error(
    Population$new(
      population = input_pop,
      treatments = character(),
      parameters = parameters,
      simulate_at_initialization = FALSE
    ),
    "non-empty list"
  )
  
  expect_error(
    Population$new(
      population = input_pop,
      treatments = list(
        c("prrt", "everolimus", NA)
      ),
      parameters = parameters,
      simulate_at_initialization = FALSE
    ),
    "uniquely named"
  )
  
  expect_error(
    Population$new(
      population = input_pop,
      treatments = list(
        bad = c("prrt", "unknown_drug", NA)
      ),
      parameters = parameters,
      simulate_at_initialization = FALSE
    ),
    "Invalid treatment"
  )
  
  expect_error(
    Population$new(
      population = input_pop,
      treatments = list(
        bad = c("prrt", NA, "everolimus")
      ),
      parameters = parameters,
      simulate_at_initialization = FALSE
    ),
    "only at the end"
  )
  
  expect_error(
    Population$new(
      population = input_pop,
      treatments = list(
        bad = c("prrt", "prrt", NA)
      ),
      parameters = parameters,
      simulate_at_initialization = FALSE
    ),
    "cannot occur more than once"
  )
})

test_that("All-cause death equals earliest valid cause-specific death", {
  P <- Population$new(
    population = input_pop,
    treatments = multiple_treatment_sequences,
    parameters = parameters,
    simulate_at_initialization = TRUE
  )
  
  simulated <- P$getPopulation()
  
  expected <- simulated[
    ,
    pmin(
      age_dead_from_other_causes,
      age_dead_from_disease_progression,
      age_dead_from_mds_aml,
      na.rm = TRUE
    )
  ]
  
  expect_equal(
    simulated$age_dead_from_all_causes,
    expected,
    tolerance = 1e-12
  )
})

test_that("All-cause death equals the earliest valid cause-specific death", {
  P <- Population$new(
    population = input_pop,
    treatments = multiple_treatment_sequences,
    parameters = parameters,
    simulate_at_initialization = TRUE
  )
  
  simulated <- P$getPopulation()
  
  expected <- simulated[
    ,
    pmin(
      age_dead_from_other_causes,
      age_dead_from_disease_progression,
      age_dead_from_mds_aml,
      na.rm = TRUE
    )
  ]
  
  expect_equal(
    simulated$age_dead_from_all_causes,
    expected,
    tolerance = 1e-12
  )
})

test_that("getOS works", {
  P1 <- Population$new(
    population = input_pop,
    treatments = multiple_treatment_sequences,
    parameters = parameters,
    simulate_at_initialization = TRUE
  )
  
  P2 <- Population$new(
    population = input_pop,
    treatments = single_treatment_sequence,
    parameters = parameters,
    simulate_at_initialization = TRUE
  )
  
  expect_no_error(P1$getOS())
  expect_no_error(P2$getOS())
  expect_identical(
    colnames(P1$getOS()[["surv_time"]]),
    colnames(P2$getOS()[["surv_time"]])
  )
  expect_identical(
    colnames(P1$getOS()[["surv_prop"]]),
    colnames(P2$getOS()[["surv_prop"]])
  )
})

test_that("getPFS works", {
  P1 <- Population$new(
    population = input_pop,
    treatments = multiple_treatment_sequences,
    parameters = parameters,
    simulate_at_initialization = TRUE
  )
  
  P2 <- Population$new(
    population = input_pop,
    treatments = single_treatment_sequence,
    parameters = parameters,
    simulate_at_initialization = TRUE
  )
  
  expect_no_error(P1$getPFS())
  expect_no_error(P2$getPFS())
  expect_identical(
    colnames(P1$getPFS()[["surv_time"]]),
    colnames(P2$getPFS()[["surv_time"]])
    )
  expect_identical(
    colnames(P1$getPFS()[["surv_prop"]]),
    colnames(P2$getPFS()[["surv_prop"]])
  )
})
  
test_that("getCumIncidenceDeathCauses works", {
  P1 <- Population$new(
    population = input_pop,
    treatments = multiple_treatment_sequences,
    parameters = parameters,
    simulate_at_initialization = TRUE
  )
  
  P2 <- Population$new(
    population = input_pop,
    treatments = single_treatment_sequence,
    parameters = parameters,
    simulate_at_initialization = TRUE
  )
  
  expect_no_error(P1$getCumIncidenceDeathCauses())
  expect_no_error(P2$getCumIncidenceDeathCauses())
  expect_identical(
    colnames(P1$getCumIncidenceDeathCauses()),
    colnames(P2$getCumIncidenceDeathCauses())
  )
})

test_that("getCumIncidenceMDSAML works", {
  P1 <- Population$new(
    population = input_pop,
    treatments = multiple_treatment_sequences,
    parameters = parameters,
    simulate_at_initialization = TRUE
  )
  
  P2 <- Population$new(
    population = input_pop,
    treatments = single_treatment_sequence,
    parameters = parameters,
    simulate_at_initialization = TRUE
  )
  
  expect_no_error(P1$getCumIncidenceMDSAML())
  expect_no_error(P2$getCumIncidenceMDSAML())
  expect_identical(
    colnames(P1$getCumIncidenceMDSAML()),
    colnames(P2$getCumIncidenceMDSAML())
  )
})

test_that("getTransitionTimes works", {
  P1 <- Population$new(
    population = input_pop,
    treatments = multiple_treatment_sequences,
    parameters = parameters,
    simulate_at_initialization = TRUE
  )
  
  P2 <- Population$new(
    population = input_pop,
    treatments = single_treatment_sequence,
    parameters = parameters,
    simulate_at_initialization = TRUE
  )
  
  expect_no_error(P1$getTransitionTimes())
  expect_no_error(P2$getTransitionTimes())
  expect_identical(
    colnames(P1$getTransitionTimes()),
    colnames(P2$getTransitionTimes())
  )
})

