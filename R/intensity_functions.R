#' Intensity for progression from treatment 1 to treatment 2
#'
#' @param t Numeric vector of ages or simulation times.
#' @param lambda_everolimus_12 Everolimus transition rate.
#' @param lambda_prrt_12 PRRT transition rate.
#' @param lambda_cabozantinib_12 Cabozantinib transition rate.
#' @param prrt_start_age PRRT start age.
#' @param prrt_stop_age PRRT stop age.
#' @param everolimus_start_age Everolimus start age.
#' @param everolimus_stop_age Everolimus stop age.
#' @param cabozantinib_start_age Cabozantinib start age.
#' @param cabozantinib_stop_age Cabozantinib stop age.
#'
#' @return A numeric vector of transition intensities.
#' 
#' @export
intensity_12 <- function(t,
                         lambda_everolimus_12,
                         lambda_prrt_12,
                         lambda_cabozantinib_12,
                         prrt_start_age,
                         prrt_stop_age,
                         everolimus_start_age,
                         everolimus_stop_age,
                         cabozantinib_start_age,
                         cabozantinib_stop_age) {
  ((t >= prrt_start_age) & (t < prrt_stop_age)) * lambda_prrt_12 +
    ((t >= everolimus_start_age) & (t < everolimus_stop_age)) * lambda_everolimus_12 +
    ((t >= cabozantinib_start_age) & (t < cabozantinib_stop_age)) * lambda_cabozantinib_12
}

#' Intensity for progression from treatment 2 to treatment 3
#'
#' @param t Numeric vector of ages or simulation times.
#' @param lambda_everolimus_23 Everolimus transition rate.
#' @param lambda_prrt_23 PRRT transition rate.
#' @param lambda_cabozantinib_23 Cabozantinib transition rate.
#' @param prrt_start_age PRRT start age.
#' @param prrt_stop_age PRRT stop age.
#' @param everolimus_start_age Everolimus start age.
#' @param everolimus_stop_age Everolimus stop age.
#' @param cabozantinib_start_age Cabozantinib start age.
#' @param cabozantinib_stop_age Cabozantinib stop age.
#'
#' @return A numeric vector of transition intensities.
#'  
#' @export
intensity_23 <- function(t,
                         lambda_everolimus_23,
                         lambda_prrt_23,
                         lambda_cabozantinib_23,
                         prrt_start_age,
                         prrt_stop_age,
                         everolimus_start_age,
                         everolimus_stop_age,
                         cabozantinib_start_age,
                         cabozantinib_stop_age) {
  ((t >= prrt_start_age) & (t < prrt_stop_age)) * lambda_prrt_23 +
    ((t >= everolimus_start_age) & (t < everolimus_stop_age)) * lambda_everolimus_23 +
    ((t >= cabozantinib_start_age) & (t < cabozantinib_stop_age)) * lambda_cabozantinib_23
}

#' Intensity for progression from treatment 3 to salvage therapy
#'
#'
#' @param t Numeric vector of ages or simulation times.
#' @param lambda_everolimus_34 Everolimus transition rate.
#' @param lambda_prrt_34 PRRT transition rate.
#' @param lambda_cabozantinib_34 Cabozantinib transition rate.
#' @param prrt_start_age PRRT start age.
#' @param prrt_stop_age PRRT stop age.
#' @param everolimus_start_age Everolimus start age.
#' @param everolimus_stop_age Everolimus stop age.
#' @param cabozantinib_start_age Cabozantinib start age.
#' @param cabozantinib_stop_age Cabozantinib stop age.
#'
#' @return A numeric vector of transition intensities.
#' 
#' @export
intensity_34 <- function(t,
                         lambda_everolimus_34,
                         lambda_prrt_34,
                         lambda_cabozantinib_34,
                         prrt_start_age,
                         prrt_stop_age,
                         everolimus_start_age,
                         everolimus_stop_age,
                         cabozantinib_start_age,
                         cabozantinib_stop_age) {
  ((t >= prrt_start_age) & (t < prrt_stop_age)) * lambda_prrt_34 +
    ((t >= everolimus_start_age) & (t < everolimus_stop_age)) * lambda_everolimus_34 +
    ((t >= cabozantinib_start_age) & (t < cabozantinib_stop_age)) * lambda_cabozantinib_34
}

#' Intensity for death following disease progression
#'
#' @param t Numeric vector of ages or simulation times.
#' @param lambda_46 lambda for death from progression
#'
#' @return A numeric vector of transition intensities.
#'  
#' @export
intensity_46 <- function(t, lambda_46) {
  t * 0 + lambda_46
}

#' Intensity for death following MDS/AML
#'
#' @param t Numeric vector of ages or simulation times.
#' @param lambda_57 lambda for death from MDS/AML
#'
#' @return A numeric vector of transition intensities.
#'  
#' @export
intensity_57 <- function(t, lambda_57) {
  t * 0 + lambda_57
}

#' Intensity for treatment-related MDS/AML
#'
#' @description
#' The intensity is proportional to the number of PRRT doses received. Doses
#' are given every `interdose_interval` years from `prrt_start_age`, up to
#' `max_doses`. No further doses are given after `prrt_stop_age`, so a person
#' who stops PRRT early keeps the risk of the doses received until then.
#'
#' @param t Numeric vector of ages or simulation times.
#' @param lambda_step lambda for progression from PRRT initiation to MDS/AML
#' @param prrt_start_age start age for prrt
#' @param prrt_stop_age stop age for prrt
#' @param max_doses maximum number of doses
#' @param interdose_interval interval between doses in years
#' 
#' @return A numeric vector of transition intensities.
#'  
#' @export
intensity_dot5 <- function(
    t, 
    lambda_step, 
    prrt_start_age, 
    prrt_stop_age = Inf,
    max_doses = 4, 
    interdose_interval = 1 / 6
    ) {
  dt_on_prrt <- pmax(pmin(t, prrt_stop_age) - prrt_start_age, 0)
  doses_received <- pmin(ceiling(dt_on_prrt / interdose_interval), max_doses)
  return(lambda_step * doses_received)
}
