#' Simulated NEC population
#'
#' An R6 class for simulating counterfactual treatment sequences, disease
#' progression, survival, and treatment-related MDS/AML outcomes.
#'
#' Each person is simulated under one or more treatment sequences. The
#' resulting population is stored in long format, with one row per person and
#' treatment sequence.
#'
#' @section Public methods:
#' \describe{
#'   \item{\code{new(population, treatments, parameters,
#'   simulate_at_initialization = TRUE)}}{
#'   Creates a new \code{Population} object and optionally runs the simulation.
#'   }
#'
#'   \item{\code{update()}}{
#'   Runs mortality, treatment-transition, disease-progression, MDS/AML, and
#'   cause-specific death simulations.
#'   }
#'
#'   \item{\code{reset()}}{
#'   Resets the object so the population can be simulated again.
#'   }
#'
#'   \item{\code{getPopulation(copy = TRUE)}}{
#'   Returns the simulated population in long format.
#'   }
#'
#'   \item{\code{getTransitionTimes()}}{
#'   Returns transition-time summaries by treatment sequence.
#'   }
#'
#'   \item{\code{getOS(time_points = c(1, 5, 10))}}{
#'   Returns overall-survival summaries by treatment sequence.
#'   }
#'
#'   \item{\code{getPFS(time_points = c(1, 5, 10),
#'   MDSAMLcensored = FALSE)}}{
#'   Returns progression-free-survival summaries by treatment sequence.
#'   }
#'
#'   \item{\code{getCumIncidenceDeathCauses(time_points = c(1, 5, 10))}}{
#'   Returns cause-specific cumulative-incidence estimates.
#'   }
#'
#'   \item{\code{getCumIncidenceMDSAML(time_points = c(1, 5, 10))}}{
#'   Returns MDS/AML cumulative-incidence estimates.
#'   }
#' }
Population <- R6::R6Class(
  classname = "Population",
  public = list(
    #' Population object
    #'
    #' @description
    #' Initializes a population, validates treatment sequences and parameters,
    #' and optionally runs the simulation.
    #'
    #' @param population A `data.table` containing at least `id`, `sex`,
    #'   `birth_cohort`, `spawn_age`, and `max_simulation_age`.
    #' @param treatments A uniquely named list of treatment sequences. Each
    #'   sequence must contain two or three treatments selected from `"prrt"`,
    #'   `"everolimus"`, and `"cabozantinib"`. Missing third treatments should
    #'   be represented by `NA`.
    #' @param simulate_at_initialization Logical. If `TRUE`, run the simulation
    #'   during initialization.
    #' @param parameters A named list of person-level transition-rate vectors.
    #'
    #' @return A new Population object.
    initialize = function(population,
                          treatments = treatment_sequences(),
                          simulate_at_initialization = TRUE,
                          parameters = NULL) {
      private$.setPopulation(population = population, copy = TRUE)
      private$.input_population <- copy(private$.population)
      private$.setParameters(parameters)
      private$.treatments <- private$.validateTreatments(treatments)
      private$.simulation_erased <- TRUE
      if (simulate_at_initialization) {
        self$update()
      }
      invisible(self)
    },
    
    #' Run the population simulation
    #'
    #' @description
    #' Generates mortality, treatment transitions, disease progression,
    #' MDS/AML onset, and cause-specific death outcomes.
    #'
    #' @return The Population object, invisibly.
    update = function() {
      stopifnot(private$.simulation_erased)
      # Generated once per person
      private$.updateParameters()
      private$.generate_death_from_other_causes()
      
      # Expand to one row per person and treatment sequence.
      private$.assignTreatments()
      
      # Treatment-dependent transitions are generated on the long data.
      private$.updateProgression12()
      private$.updateProgression23()
      private$.updateProgression34()
      private$.updateProgression46()
      private$.updateProgressionDot5()
      private$.updateProgression57()
      
      private$.simulation_erased <- FALSE
      
      invisible(self)
    },
    #' Reset the simulation
    #'
    #' @description
    #' Removes simulated outcomes and marks the object as ready to simulate
    #' again.
    #'
    #' @return The Population object, invisibly.
    reset = function() {
      # Restore the one-row-per-person input population, so update() can
      # simulate it again from scratch.
      private$.setPopulation(population = private$.input_population, copy = TRUE)
      private$.simulation_erased <- TRUE
      invisible(self)
    },
    #' Return the simulated population
    #'
    #' @param copy Logical. If `TRUE`, return a copy of the population.
    #'   If `FALSE`, return the internal population by reference.
    #'
    #' @return A `data.table` in long format with one row per person and
    #'   treatment sequence.
    getPopulation = function(copy = TRUE) {
      if (copy) {
        return(copy(private$.population))
      } else {
        return(private$.population)
      }
    },
    #' Summarize transition times
    #'
    #' @description
    #' Summarizes transition frequencies and durations by treatment sequence.
    #'
    #' @return A `data.table` with transition summaries by treatment sequence.
    getTransitionTimes = function() {
      private$.summariseTransitionTimes(dt = private$.population)
    },
    #' Estimate overall survival
    #'
    #' @param time_points Numeric vector of follow-up times, in years, at which
    #'   survival probabilities are reported.
    #'
    #' @return A list containing `surv_prop` and `surv_time`, stratified by
    #'   treatment sequence.
    getOS = function(time_points = c(1, 5, 10)) {
      private$.population[, `:=`(event = 1, time_to_event = age_dead_from_all_causes - spawn_age)]
      
      OS <- private$.getKM(dt = private$.population, time_points = time_points)
      
      private$.population[, c("event", "time_to_event") := NULL]
      
      return(OS)
    },
    #' Estimate progression-free survival
    #'
    #' @param time_points Numeric vector of follow-up times, in years, at which
    #'   survival probabilities are reported.
    #' @param MDSAMLcensored Logical. If `FALSE` (default), PFS is the time to
    #'   the earliest of salvage treatment initiation, MDS/AML, and death. If
    #'   `TRUE`, MDS/AML onset is treated as a censoring event instead.
    #'
    #' @return A list containing `surv_prop` and `surv_time`, stratified by
    #'   treatment sequence.
    getPFS = function(time_points = c(1, 5, 10),
                      MDSAMLcensored = FALSE
                      ) {
      if (isTRUE(MDSAMLcensored)) {
        private$.population[, event := ifelse(is.na(age_with_mds_aml), 1, 0)]
        private$.population[, time_to_event := pmin(
          age_with_mds_aml,
          age_progression_to_salvage_therapy,
          age_dead_from_all_causes,
          na.rm = TRUE
        ) -
          spawn_age]
      } else {
        private$.population[, event := as.numeric(1)]
        private$.population[, time_to_event := pmin(age_progression_to_salvage_therapy,
                                                    age_with_mds_aml,
                                                    age_dead_from_all_causes,
                                                    na.rm = TRUE) -
                              spawn_age]
      }
      
      PFS <- private$.getKM(dt = private$.population, time_points = time_points)
      
      cols_to_remove <- c("event", "time_to_event")
      private$.population[, (intersect(names(private$.population), cols_to_remove)) := NULL]
      return(PFS)
    },
    #' Estimate cumulative incidence of causes of death
    #'
    #' @param time_points Numeric vector of follow-up times, in years.
    #'
    #' @return A `data.table` containing cumulative incidence estimates by
    #'   treatment sequence, death cause, and time.
    getCumIncidenceDeathCauses = function(time_points = c(1, 5, 10)) {
      # browser()
      death_cause_labels <- c(
        "age_dead_from_other_causes" = "Other causes",
        "age_dead_from_mds_aml" = "MDS/AML",
        "age_dead_from_disease_progression" = "Disease progression"
      )
      # Earliest cause-specific death age; missing causes never win.
      death_ages <- as.matrix(private$.population[, names(death_cause_labels), with = FALSE])
      death_ages[is.na(death_ages)] <- Inf
      private$.population[, `:=`(event = as.numeric(max.col(-death_ages, ties.method = "first")),
                                 time_to_event = age_dead_from_all_causes - spawn_age)]
      
      cif_death <- private$.getCumIncidence(dt = private$.population, time_points = time_points)
      
      cif_death[, event := factor(event, levels = 1:3, labels = death_cause_labels)]
      
      private$.population[, `:=`(event = NULL, time_to_event = NULL)]
      
      return(cif_death)
    },
    #' Estimate cumulative incidence of MDS/AML
    #'
    #' @param time_points Numeric vector of follow-up times, in years.
    #'
    #' @return A `data.table` containing MDS/AML cumulative incidence estimates
    #'   by treatment sequence and time.
    getCumIncidenceMDSAML = function(time_points = c(1, 5, 10)) {
      # browser()
      # Death before MDS/AML is a competing event (2), not censoring.
      private$.population[, `:=`(
        event = fcase(
          !is.na(age_with_mds_aml), 1,
          !is.na(age_dead_from_all_causes), 2,
          default = 0
        ),
        time_to_event = pmin(age_with_mds_aml, age_dead_from_all_causes, na.rm = TRUE) -
          spawn_age
      )]

      cif_mdsaml <- private$.getCumIncidence(dt = private$.population, time_points = time_points)

      cif_mdsaml <- cif_mdsaml[event == "1"]
      cif_mdsaml[, event := "MDS/AML"]
      
      private$.population[, `:=`(event = NULL, time_to_event = NULL)]
      
      return(cif_mdsaml)
    }
  ),
  private = list(
    .population = NULL,
    .input_population = NULL,
    .simulation_erased = NULL,
    .treatments = NULL,
    .parameters = NULL,
    .setPopulation = function(population, copy) {
      stopifnot(is.data.table(population))
      stopifnot(all(
        c(
          "id",
          "sex",
          "birth_cohort",
          "spawn_age",
          "max_simulation_age"
        ) %in%
          names(population)
      ))
      if (copy) {
        private$.population <- copy(population)
      } else {
        private$.population <- population
      }
      setkey(private$.population, id)
      invisible(self)
    },
    .setParameters = function(parameters) {
      stopifnot(!is.null(parameters))
      stopifnot(is.list(parameters))
      stopifnot(all(
        c(
          "lambda_everolimus_12",
          "lambda_everolimus_23",
          "lambda_everolimus_34",
          "lambda_prrt_12",
          "lambda_prrt_23",
          "lambda_prrt_34",
          "lambda_cabozantinib_12",
          "lambda_cabozantinib_23",
          "lambda_cabozantinib_34",
          "lambda_46",
          "lambda_57",
          "lambda_dot5"
        ) %in%
          names(parameters)
      ))
      private$.parameters <- parameters
    },
    .updateParameters = function() {
      # browser()
      private$.population[, `:=`(
        params_lambda_everolimus_12 = private$.parameters$lambda_everolimus_12,
        params_lambda_everolimus_23 = private$.parameters$lambda_everolimus_23,
        params_lambda_everolimus_34 = private$.parameters$lambda_everolimus_34,
        params_lambda_prrt_12 = private$.parameters$lambda_prrt_12,
        params_lambda_prrt_23 = private$.parameters$lambda_prrt_23,
        params_lambda_prrt_34 = private$.parameters$lambda_prrt_34,
        params_lambda_cabozantinib_12 = private$.parameters$lambda_cabozantinib_12,
        params_lambda_cabozantinib_23 = private$.parameters$lambda_cabozantinib_23,
        params_lambda_cabozantinib_34 = private$.parameters$lambda_cabozantinib_34,
        params_lambda_46 = private$.parameters$lambda_46,
        params_lambda_57 = private$.parameters$lambda_57,
        params_lambda_dot5 = private$.parameters$lambda_dot5
      )]
      invisible(self)
    },
    # Transition rate of the treatment given in `treatment_var`, one per row.
    # Used as the thinning majorizer: with a single active treatment the
    # intensity is constant and every candidate point is accepted.
    .activeTreatmentRate = function(dt, treatment_var, transition) {
      tx <- dt[[treatment_var]]
      rate <- fcase(
        tx == "prrt", dt[[paste0("params_lambda_prrt_", transition)]],
        tx == "everolimus", dt[[paste0("params_lambda_everolimus_", transition)]],
        tx == "cabozantinib", dt[[paste0("params_lambda_cabozantinib_", transition)]]
      )
      stopifnot(!anyNA(rate))
      rate
    },
    .validateTreatments = function(treatments) {
      # browser()
      allowed_treatments <- c("prrt", "everolimus", "cabozantinib")
      
      if (!is.list(treatments) || length(treatments) == 0L) {
        stop("`treatments` must be a non-empty list.")
      }
      
      if (is.null(names(treatments)) ||
          any(names(treatments) == "") ||
          anyDuplicated(names(treatments))) {
        stop("`treatments` must be a uniquely named list. ",
             "The names will be stored in `tx_seq`.")
      }
      
      treatments <- lapply(treatments, function(x) {
        x <- as.character(x)
        if (length(x) > 3L) {
          stop("Each treatment sequence may contain at most three treatments.")
        }
        
        # Extend shorter sequences to three positions.
        length(x) <- 3L
        x[is.na(x) | x == ""] <- NA_character_
        
        invalid <- setdiff(stats::na.omit(x), allowed_treatments)
        
        if (length(invalid) > 0L) {
          stop("Invalid treatment(s): ", paste(invalid, collapse = ", "))
        }
        
        if (anyDuplicated(stats::na.omit(x))) {
          stop("A treatment cannot occur more than once in the same sequence.")
        }
        
        # Once an NA appears, later lines must also be NA.
        nonmissing_positions <- which(!is.na(x))
        
        if (length(nonmissing_positions) > 0L &&
            !identical(nonmissing_positions, seq_len(max(nonmissing_positions)))) {
          stop("Missing treatments are allowed only at the end of a sequence.")
        }
        x
      })
      treatments
    },
    .assignTreatments = function() {
      # browser()
      sequence_table <- rbindlist(lapply(seq_along(private$.treatments), function(i) {
        sequence <- private$.treatments[[i]]
        
        data.table(
          tx_seq_id = i,
          tx_seq = names(private$.treatments)[i],
          treatment1 = sequence[1],
          treatment2 = sequence[2],
          treatment3 = sequence[3]
        )
      }))
      
      # Cartesian join: each person receives every treatment sequence.
      private$.population[, join_id_ := 1L]
      sequence_table[, join_id_ := 1L]
      
      private$.population <- merge(
        private$.population,
        sequence_table,
        by = "join_id_",
        allow.cartesian = TRUE,
        sort = FALSE
      )
      
      private$.population[, join_id_ := NULL]
      
      setorder(private$.population, tx_seq_id, id)
      
      private$.population[, `:=`(
        prrt_start_age = Inf,
        prrt_stop_age = Inf,
        everolimus_start_age = Inf,
        everolimus_stop_age = Inf,
        cabozantinib_start_age = Inf,
        cabozantinib_stop_age = Inf
      )]
      
      treatment_names <- c("prrt", "everolimus", "cabozantinib")
      
      # First-line treatment starts at model entry.
      for (tx in treatment_names) {
        start_var <- paste0(tx, "_start_age")
        stop_var <- paste0(tx, "_stop_age")
        
        private$.population[treatment1 == tx, (start_var) := spawn_age]
        
        private$.population[treatment1 == tx, (stop_var) := pmin(age_dead_from_other_causes, max_simulation_age, na.rm = TRUE)]
      }
      
      setkey(private$.population, tx_seq_id, id)
      
      invisible(self)
    },
    .generate_death_from_other_causes = function() {
      # browser()
      setorder(private$.population, id)
      private$.population[, cohort := birth_cohort]
      private$.population[, age_dead_from_other_causes := mortality::vdraw_age_dead_other_causes_in_years(
        start_age_in_years = spawn_age,
        stop_age_in_years = max_simulation_age,
        population = private$.population,
        region = "USA",
        source = "HMD"
      )]
      private$.population[, cohort := NULL]
      invisible(self)
    },
    .updateTreatmentAges = function(transition_age_var,
                                    stop_treatment_var = NULL,
                                    start_treatment_var = NULL) {
      age <- private$.population[[transition_age_var]]
      
      # Stop current treatment.
      if (!is.null(stop_treatment_var)) {
        stop_tx <- private$.population[[stop_treatment_var]]
        
        private$.population[, `:=`(
          prrt_stop_age = fifelse(
            !is.na(age) & !is.na(stop_tx) & stop_tx == "prrt",
            age,
            prrt_stop_age
          ),
          everolimus_stop_age = fifelse(
            !is.na(age) & !is.na(stop_tx) & stop_tx == "everolimus",
            age,
            everolimus_stop_age
          ),
          cabozantinib_stop_age = fifelse(
            !is.na(age) & !is.na(stop_tx) & stop_tx == "cabozantinib",
            age,
            cabozantinib_stop_age
          )
        )]
      }
      
      # Start next treatment and provisionally continue it until other-cause death.
      if (!is.null(start_treatment_var)) {
        start_tx <- private$.population[[start_treatment_var]]
        
        private$.population[, `:=`(
          prrt_start_age = fifelse(
            !is.na(age) & !is.na(start_tx) & start_tx == "prrt",
            age,
            prrt_start_age
          ),
          everolimus_start_age = fifelse(
            !is.na(age) & !is.na(start_tx) & start_tx == "everolimus",
            age,
            everolimus_start_age
          ),
          cabozantinib_start_age = fifelse(
            !is.na(age) & !is.na(start_tx) & start_tx == "cabozantinib",
            age,
            cabozantinib_start_age
          ),
          prrt_stop_age = fifelse(
            !is.na(age) & !is.na(start_tx) & start_tx == "prrt",
            age_dead_from_other_causes,
            prrt_stop_age
          ),
          everolimus_stop_age = fifelse(
            !is.na(age) & !is.na(start_tx) & start_tx == "everolimus",
            age_dead_from_other_causes,
            everolimus_stop_age
          ),
          cabozantinib_stop_age = fifelse(
            !is.na(age) & !is.na(start_tx) & start_tx == "cabozantinib",
            age_dead_from_other_causes,
            cabozantinib_stop_age
          )
        )]
      }
      
      invisible(self)
    },
    # Age of progression to "2. Stable (Tx 2)"
    .updateProgression12 = function() {
      # One row per person and first treatment.
      history_1 <- unique(private$.population[, .(
        id,
        treatment1,
        spawn_age,
        max_simulation_age,
        age_dead_from_other_causes,
        
        params_lambda_everolimus_12,
        params_lambda_prrt_12,
        params_lambda_cabozantinib_12,
        
        prrt_start_age,
        prrt_stop_age,
        everolimus_start_age,
        everolimus_stop_age,
        cabozantinib_start_age,
        cabozantinib_stop_age
      )], by = c("id", "treatment1"))
      
      if (nrow(history_1) == 0L) {
        return(invisible(self))
      }
      
      # Intensity evaluated on the reduced first-line history table.
      l <- function(t, pop = history_1, ...) {
        intensity_12(
          t = t,
          lambda_everolimus_12 =
            pop$params_lambda_everolimus_12,
          lambda_prrt_12 =
            pop$params_lambda_prrt_12,
          lambda_cabozantinib_12 =
            pop$params_lambda_cabozantinib_12,
          
          prrt_start_age = pop$prrt_start_age,
          prrt_stop_age = pop$prrt_stop_age,
          
          everolimus_start_age =
            pop$everolimus_start_age,
          everolimus_stop_age =
            pop$everolimus_stop_age,
          
          cabozantinib_start_age =
            pop$cabozantinib_start_age,
          cabozantinib_stop_age =
            pop$cabozantinib_stop_age
        )
      }
      
      # Only treatment 1 is active on [t_min, t_max], so its rate is an exact majorizer.
      l_star <- matrix(private$.activeTreatmentRate(history_1, "treatment1", "12"), ncol = 1)
      
      history_1[, age_progression_to_stable_2 :=
                  nhppp::vdraw_intensity(
                    lambda = l,
                    lambda_maj_matrix = l_star,
                    t_min = spawn_age,
                    t_max = age_dead_from_other_causes,
                    rate_matrix_t_min = min(spawn_age),
                    rate_matrix_t_max = max(max_simulation_age),
                    atmost1 = TRUE
                  )]
      
      # Join the shared first-line transition back to all strategies
      # with the same person and first-line treatment.
      private$.population[history_1[, .(id, treatment1, age_progression_to_stable_2)], on = .(id, treatment1), age_progression_to_stable_2 :=
                            i.age_progression_to_stable_2]
      
      # Stop treatment 1 and start treatment 2.
      private$.updateTreatmentAges(
        transition_age_var =
          "age_progression_to_stable_2",
        stop_treatment_var = "treatment1",
        start_treatment_var = "treatment2"
      )
      
      invisible(self)
    },
    # Age of progression to "3. Stable (Tx 3)"
    .updateProgression23 = function() {
      eligible <- private$.population[!is.na(treatment2) &
                                        !is.na(age_progression_to_stable_2)]
      
      if (nrow(eligible) == 0L) {
        return(invisible(self))
      }
      
      # One row per person and complete treatment history through line 2.
      history_2 <- unique(eligible[, .(
        id,
        treatment1,
        treatment2,
        
        spawn_age,
        max_simulation_age,
        age_dead_from_other_causes,
        age_progression_to_stable_2,
        
        params_lambda_everolimus_23,
        params_lambda_prrt_23,
        params_lambda_cabozantinib_23,
        
        prrt_start_age,
        prrt_stop_age,
        everolimus_start_age,
        everolimus_stop_age,
        cabozantinib_start_age,
        cabozantinib_stop_age
      )], by = c("id", "treatment1", "treatment2"))
      
      l <- function(t, pop = history_2, ...) {
        intensity_23(
          t = t,
          lambda_everolimus_23 =
            pop$params_lambda_everolimus_23,
          lambda_prrt_23 =
            pop$params_lambda_prrt_23,
          lambda_cabozantinib_23 =
            pop$params_lambda_cabozantinib_23,
          
          prrt_start_age = pop$prrt_start_age,
          prrt_stop_age = pop$prrt_stop_age,
          
          everolimus_start_age =
            pop$everolimus_start_age,
          everolimus_stop_age =
            pop$everolimus_stop_age,
          
          cabozantinib_start_age =
            pop$cabozantinib_start_age,
          cabozantinib_stop_age =
            pop$cabozantinib_stop_age
        )
      }
      
      range_t <- history_2[, range(spawn_age, max_simulation_age, na.rm = TRUE)]
      
      l_star <- matrix(private$.activeTreatmentRate(history_2, "treatment2", "23"), ncol = 1)
      
      history_2[, age_after_treatment2 :=
                  nhppp::vdraw_intensity(
                    lambda = l,
                    lambda_maj_matrix = l_star,
                    t_min = age_progression_to_stable_2,
                    t_max = age_dead_from_other_causes,
                    rate_matrix_t_min = range_t[1],
                    rate_matrix_t_max = range_t[2],
                    atmost1 = TRUE
                  )]
      
      # Join the shared second-line event back to all strategies
      # with the same treatment history.
      private$.population[history_2[, .(id, treatment1, treatment2, age_after_treatment2)], on = .(id, treatment1, treatment2), age_after_treatment2 := i.age_after_treatment2]
      
      # Initialize the destination-state column.
      private$.population[, age_progression_to_stable_3 := NA_real_]
      
      # Three-drug strategies move to stable state 3.
      private$.population[!is.na(treatment3) &
                            !is.na(age_after_treatment2), age_progression_to_stable_3 := age_after_treatment2]
      
      # Two-drug strategies move directly to salvage.
      # Their stable-3 age stays NA.
      private$.population[is.na(treatment3) &
                            !is.na(age_after_treatment2), age_progression_to_salvage_therapy := age_after_treatment2]
      
      # Treatment 2 stops at the common transition age.
      private$.updateTreatmentAges(transition_age_var = "age_after_treatment2",
                                   stop_treatment_var = "treatment2")
      
      # Treatment 3 starts only in three-drug strategies.
      private$.updateTreatmentAges(transition_age_var = "age_progression_to_stable_3",
                                   start_treatment_var = "treatment3")
      
      # No longer needed after routing the event.
      private$.population[, age_after_treatment2 := NULL]
      
      invisible(self)
    },
    # Progression from "3. Stable (Tx 3)" (sometimes present, sometimes not) to "4. Salvage"
    .updateProgression34 = function() {
      # No third treatment: move directly to salvage.
      eligible <- private$.population[!is.na(treatment3) &
                                        !is.na(age_progression_to_stable_3)]
      
      if (nrow(eligible) == 0L) {
        return(invisible(self))
      }
      
      l <- function(t, pop = eligible, ...) {
        intensity_34(
          t = t,
          lambda_everolimus_34 = pop$params_lambda_everolimus_34,
          lambda_prrt_34 = pop$params_lambda_prrt_34,
          lambda_cabozantinib_34 = pop$params_lambda_cabozantinib_34,
          prrt_start_age = pop$prrt_start_age,
          prrt_stop_age = pop$prrt_stop_age,
          everolimus_start_age = pop$everolimus_start_age,
          everolimus_stop_age = pop$everolimus_stop_age,
          cabozantinib_start_age = pop$cabozantinib_start_age,
          cabozantinib_stop_age = pop$cabozantinib_stop_age
        )
      }
      
      range_t <- eligible[, range(spawn_age, max_simulation_age)]
      
      l_star <- matrix(private$.activeTreatmentRate(eligible, "treatment3", "34"), ncol = 1)
      
      private$.population[!is.na(treatment3) &
                            !is.na(age_progression_to_stable_3), age_progression_to_salvage_therapy := nhppp::vdraw_intensity(
                              lambda = l,
                              lambda_maj_matrix = l_star,
                              t_min = age_progression_to_stable_3,
                              t_max = age_dead_from_other_causes,
                              rate_matrix_t_min = range_t[1],
                              rate_matrix_t_max = range_t[2],
                              atmost1 = TRUE,
                              atmostB = 20
                            )]
      
      private$.updateTreatmentAges(transition_age_var = "age_progression_to_salvage_therapy",
                                   stop_treatment_var = "treatment3")
      
      invisible(self)
    },
    # Age of progression from "4. Salvage" to "6. Death (progression)"
    .updateProgression46 = function() {
      # browser()
      private$.population[, age_dead_from_disease_progression := NA_real_]

      eligible <- private$.population[!is.na(age_progression_to_salvage_therapy)]

      if (nrow(eligible) > 0L) {
        l <- function(t, pop = eligible, ...) {
          intensity_46(t = t, lambda_46 = pop$params_lambda_46)
        }
        range_t <- eligible[, range(spawn_age, max_simulation_age)]
        l_star <- matrix(eligible$params_lambda_46, ncol = 1)
        private$.population[!is.na(age_progression_to_salvage_therapy), age_dead_from_disease_progression := nhppp::vdraw_intensity(
          lambda = l,
          lambda_maj_matrix = l_star,
          t_min = age_progression_to_salvage_therapy,
          t_max = age_dead_from_other_causes,
          rate_matrix_t_min = range_t[1],
          rate_matrix_t_max = range_t[2],
          atmost1 = TRUE
        )]
      }
      # update age_dead_from_all_causes
      private$.population[, age_dead_from_all_causes := pmin(age_dead_from_disease_progression,
                                                             age_dead_from_other_causes,
                                                             na.rm = TRUE)]
      invisible(self)
    },
    # age of progression to from 1, 2, 3, or 4 to "5. MDS/AML"
    .updateProgressionDot5 = function() {
      # browser()
      # Initialize columns for all strategies.
      private$.population[, `:=`(latent_age_with_mds_aml = NA_real_,
                                 age_with_mds_aml = NA_real_)]
      
      # Under the current intensity_dot5() formulation, risk depends on
      # PRRT start and stop ages (doses stop accruing at PRRT stop).
      eligible <- private$.population[is.finite(prrt_start_age)]
      
      if (nrow(eligible) == 0L) {
        return(invisible(self))
      }
      
      # One row per person and unique PRRT start and stop ages.
      # Strategies with the same id and PRRT history share the same
      # latent MDS/AML process.
      prrt_histories <- unique(eligible[, .(id,
                                            prrt_start_age,
                                            prrt_stop_age,
                                            age_dead_from_other_causes,
                                            params_lambda_dot5)], by = c("id", "prrt_start_age", "prrt_stop_age"))
      
      # Intensity evaluated on unique PRRT histories.
      l <- function(t, pop = prrt_histories, ...) {
        intensity_dot5(
          t = t,
          lambda_step = pop$params_lambda_dot5,
          prrt_start_age = pop$prrt_start_age,
          prrt_stop_age = pop$prrt_stop_age
        )
      }
      
      range_t <- prrt_histories[, range(prrt_start_age, age_dead_from_other_causes, na.rm = TRUE)]
      
      # The intensity is non-decreasing and flat after PRRT stops, so its
      # value at the PRRT stop age is the maximum.
      l_star <- matrix(l(prrt_histories$prrt_stop_age), ncol = 1)
      stopifnot(all(is.finite(l_star)))
      
      # Simulate one latent MDS/AML age per unique PRRT history.
      prrt_histories[, latent_age_with_mds_aml :=
                       nhppp::vdraw_intensity(
                         lambda = l,
                         lambda_maj_matrix = l_star,
                         t_min = prrt_start_age,
                         t_max = age_dead_from_other_causes,
                         rate_matrix_t_min = range_t[1],
                         rate_matrix_t_max = range_t[2],
                         atmost1 = TRUE,
                         atmostB = 20
                       )]
      
      # Join the shared latent event back to all strategies with the same
      # person and PRRT history.
      private$.population[prrt_histories[, .(id, prrt_start_age, prrt_stop_age, latent_age_with_mds_aml)], on = .(id, prrt_start_age, prrt_stop_age), latent_age_with_mds_aml := i.latent_age_with_mds_aml]
      
      # The shared latent event is observed only if it occurs before
      # strategy-specific all-cause death.
      private$.population[, age_with_mds_aml := fifelse(
        !is.na(latent_age_with_mds_aml) &
          latent_age_with_mds_aml < age_dead_from_all_causes,
        latent_age_with_mds_aml,
        NA_real_
      )]
      
      # Remove transitions precluded by MDS/AML.
      competing_event_vars <- c(
        "age_progression_to_stable_2",
        "age_progression_to_stable_3",
        "age_progression_to_salvage_therapy",
        "age_dead_from_disease_progression"
      )
      
      for (var in competing_event_vars) {
        private$.population[!is.na(age_with_mds_aml) &
                              !is.na(get(var)) &
                              get(var) >= age_with_mds_aml, (var) := NA_real_]
      }
      
      treatment_names <- c("prrt", "everolimus", "cabozantinib")
      treatment_start_age_vars <- paste0(treatment_names, "_start_age")
      treatment_stop_age_vars <- paste0(treatment_names, "_stop_age")
      
      # Any treatment scheduled to start at or after MDS/AML does not start.
      for (start_var in treatment_start_age_vars) {
        private$.population[!is.na(age_with_mds_aml) &
                              get(start_var) >= age_with_mds_aml, (start_var) := Inf]
      }
      
      # Stop an active treatment at MDS/AML. If a treatment never started,
      # its stop age is also reset to Inf.
      for (i in seq_along(treatment_names)) {
        start_var <- treatment_start_age_vars[i]
        stop_var <- treatment_stop_age_vars[i]
        
        private$.population[!is.na(age_with_mds_aml) &
                              is.finite(get(start_var)) &
                              get(start_var) < age_with_mds_aml &
                              get(stop_var) >= age_with_mds_aml, (stop_var) := age_with_mds_aml]
        
        private$.population[!is.na(age_with_mds_aml) &
                              is.infinite(get(start_var)), (stop_var) := Inf]
      }
      
      invisible(self)
    },
    # Age of progression from "5. MDS/AML" to "6. Death (MDS/AML)"
    .updateProgression57 = function() {
      # Initialize for everyone.
      private$.population[, age_dead_from_mds_aml := NA_real_]
      
      eligible <- private$.population[!is.na(age_with_mds_aml)]
      
      if (nrow(eligible) == 0L) {
        private$.population[, age_dead_from_all_causes := pmin(age_dead_from_other_causes,
                                                               age_dead_from_disease_progression,
                                                               na.rm = TRUE)]
        
        return(invisible(self))
      }
      
      # One row per person and shared MDS/AML onset age.
      mdsaml_histories <- unique(eligible[, .(id,
                                              age_with_mds_aml,
                                              age_dead_from_other_causes,
                                              params_lambda_57)], by = c("id", "age_with_mds_aml"))
      
      l <- function(t, pop = mdsaml_histories, ...) {
        intensity_57(t = t, lambda_57 = pop$params_lambda_57)
      }
      
      range_t <- mdsaml_histories[, range(age_with_mds_aml, age_dead_from_other_causes, na.rm = TRUE)]
      
      l_star <- matrix(mdsaml_histories$params_lambda_57, ncol = 1)
      
      # Once per shared MDS/AML history.
      mdsaml_histories[, age_dead_from_mds_aml :=
                         nhppp::vdraw_intensity(
                           lambda = l,
                           lambda_maj_matrix = l_star,
                           t_min = age_with_mds_aml,
                           t_max = age_dead_from_other_causes,
                           rate_matrix_t_min = range_t[1],
                           rate_matrix_t_max = range_t[2],
                           atmost1 = TRUE,
                           atmostB = 20
                         )]
      
      # Join shared MDS/AML death age back to all matching strategies.
      private$.population[mdsaml_histories[, .(id, age_with_mds_aml, age_dead_from_mds_aml)], on = .(id, age_with_mds_aml), age_dead_from_mds_aml :=
                            i.age_dead_from_mds_aml]
      
      # Recalculate final all-cause death from valid cause-specific deaths.
      private$.population[, age_dead_from_all_causes := pmin(
        age_dead_from_other_causes,
        age_dead_from_disease_progression,
        age_dead_from_mds_aml,
        na.rm = TRUE
      )]
      
      invisible(self)
    },
    .summariseTransitionTimes = function(dt) {
      # browser()
      vars_to_melt <- c(
        "spawn_age",
        "age_progression_to_stable_2",
        "age_progression_to_stable_3",
        "age_progression_to_salvage_therapy",
        "age_dead_from_disease_progression",
        "age_with_mds_aml",
        "age_dead_from_mds_aml",
        "age_dead_from_other_causes"
      )
      id_vars <- setdiff(names(dt), vars_to_melt)
      dt_long <- melt(
        dt,
        id.vars = id_vars,
        measure.vars = vars_to_melt,
        variable.name = "state",
        value.name = "age",
        na.rm = TRUE
      )
      setorderv(dt_long, cols = c("tx_seq", "id", "age"))
      
      dt_long[, `:=`(
        transition = paste0(state, " -> ", shift(state, n = 1, type = "lead")),
        transition_years = shift(age, n = 1, type = "lead") - age
      ), by = .(tx_seq, id)]
      dt_long <- dt_long[!grepl("dead", state), ]
      
      n_population_units <- uniqueN(dt_long, by = c("tx_seq", "id"))
      
      dt_summary <- dt_long[, .(
        N = .N,
        proportion = .N / n_population_units,
        mean = mean(transition_years, na.rm = TRUE),
        sd = sd(transition_years, na.rm = TRUE),
        median = median(transition_years, na.rm = TRUE),
        q25 = quantile(transition_years, 0.25, na.rm = TRUE),
        q75 = quantile(transition_years, 0.75, na.rm = TRUE),
        min = min(transition_years, na.rm = TRUE),
        max = max(transition_years, na.rm = TRUE)
      ), by = .(transition, tx_seq)]
      
      return(dt_summary)
    },
    .getKM = function(dt, time_points) {
      # browser()
      km_fit <- survival::survfit(survival::Surv(time = time_to_event, event = event) ~ tx_seq,
                                  data = dt)
      km_summary <- summary(km_fit, times = time_points)
      
      if (is.vector(km_summary$table)) {
        surv_time <- as.data.table(matrix(
          km_summary$table,
          nrow = length(private$.treatments),
          byrow = TRUE,
          dimnames = list(NULL, names(km_summary$table))
        ))
        surv_time[, tx := names(private$.treatments)]
        setcolorder(surv_time, "tx")
      } else {
        surv_time <- data.table(km_summary$table, keep.rownames = "tx")
        surv_time[, tx := sub("^tx_seq=", "", tx)]
      }
      
      if (is.null(km_summary$strata)) {
        km_summary$strata <- names(private$.treatments)
      } else {
        km_summary$strata <- sub("^tx_seq=", "", km_summary$strata)
      }
      
      surv_prop <- data.table(
        tx = km_summary$strata,
        time = km_summary$time,
        n.risk = as.integer(km_summary$n.risk),
        n.event = as.integer(km_summary$n.event),
        # Unrounded: these estimates feed uncertainty summaries
        survival = km_summary$surv,
        lower_ci = km_summary$lower,
        upper_ci = km_summary$upper
      )
      return(list("surv_prop" = surv_prop, "surv_time" = surv_time))
    },
    .getCumIncidence = function(dt, time_points) {
      # browser()
      cif_fit <- cmprsk::cuminc(
        ftime = dt$time_to_event,
        fstatus = dt$event,
        group = dt$tx_seq
      )
      
      cif_summary <- cmprsk::timepoints(cif_fit, times = time_points)
      
      estimates <- cif_summary$est
      variances <- cif_summary$var
      
      if (is.null(estimates) || nrow(estimates) == 0L) {
        return(
          data.table(
            tx = character(),
            event = character(),
            time = numeric(),
            n.risk = integer(),
            incidence = numeric(),
            se = numeric(),
            lower_ci = numeric(),
            upper_ci = numeric()
          )
        )
      }
      
      curve_names <- rownames(estimates)
      
      # Last token is the event code.
      curve_event <- sub("^.*\\s+([^\\s]+)$", "\\1", curve_names)
      
      curve_tx <- sub("\\s+\\d+$", "", curve_names)
      
      standard_error <- sqrt(variances)
      z <- stats::qnorm(0.975)
      lower_ci <- pmax(estimates - z * standard_error, 0)
      upper_ci <- pmin(estimates + z * standard_error, 1)
      
      # One pass over the data per time point, counted by treatment sequence.
      tx_seq <- as.character(dt$tx_seq)
      at_risk <- vapply(time_points, function(t) {
        n_at_risk <- tapply(dt$time_to_event >= t, tx_seq, sum, na.rm = TRUE)
        as.integer(n_at_risk[curve_tx])
      }, integer(length(curve_tx)))

      # Unrounded: these estimates feed sim_data() and the GAM emulator.
      data.table(
        tx = rep(curve_tx, times = ncol(estimates)),
        event = rep(curve_event, times = ncol(estimates)),
        time = rep(as.numeric(colnames(estimates)), each = nrow(estimates)),
        n.risk = as.vector(at_risk),
        incidence = as.vector(estimates),
        se = as.vector(standard_error),
        lower_ci = as.vector(lower_ci),
        upper_ci = as.vector(upper_ci)
      )
    }
  )
)
