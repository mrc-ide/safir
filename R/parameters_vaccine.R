# --------------------------------------------------------------------------------
#   infection process for vaccination model (multiple doses, no types)
#   Sean L. Wu (slwood89@gmail.com)
#   July 2021
# --------------------------------------------------------------------------------


#' @title Get vaccine efficacy and Ab titre parameters
#' @param vaccine which vaccine? should be one of: "Pfizer", "AstraZeneca", "Sinovac", "Moderna"
#' @param max_dose maximum number of doses
#' @param correlated are doses correlated?
#' @description Get parameters for vaccine efficacy and antibody titre decay rate.
#' @export
get_vaccine_ab_titre_parameters <- function(vaccine, max_dose = 2, correlated = FALSE) {
  stopifnot(max_dose %in% c(2,3))
  stopifnot(is.logical(correlated))
  stopifnot(vaccine %in% c("Pfizer", "AstraZeneca", "Sinovac", "Moderna"))

  mu_ab_list <- data.frame(name = c("Pfizer", "AstraZeneca", "Sinovac", "Moderna"),
                           mu_ab_d1 = c(13/94, 1/59, 28/164, ((185+273)/2)/321),
                           mu_ab_d2 = c(223/94, 32/59, 28/164, 654/158))

  hl_s <- 108 # Half life of antibody decay - short
  hl_l <- 3650 # Half life of antibody decay - long
  period_s <- 250
  t_period_l <- 365 # Time point at which to have switched to longest half-life
  time_to_decay <- t_period_l - period_s # time in days to reach longest half-life
  dr_s <- -log(2)/hl_s # Corresponding decay rate in days for half life above
  dr_l <- -log(2)/hl_l

  mu_ab_d1 <- mu_ab_list[mu_ab_list$name == vaccine, "mu_ab_d1"] # mean titre dose 1
  mu_ab_d2 <- mu_ab_list[mu_ab_list$name == vaccine, "mu_ab_d2"] # mean titre dose 2
  ab_50 <- 0.2 # titre relative to convalescent required to provide 50% protection from infection, on linear scale
  ab_50_severe <- 0.03
  std10 <- 0.44 # Pooled standard deviation of antibody level on log10 data
  k <- 2.94 # shape parameter of efficacy curve

  dr_vec <- c(
    rep(dr_s, period_s),
    seq(dr_s, dr_l, length.out = time_to_decay),
    dr_l
  )

  if (max_dose == 3) {
    mu_ab = rep(c(mu_ab_d1, mu_ab_d2), times = c(1,2))
  } else {
    mu_ab = rep(c(mu_ab_d1, mu_ab_d2), times = c(1,1))
  }

  parameters <- list(
    dr_vec = dr_vec,
    mu_ab = mu_ab,
    std10 = std10,
    ab_50 = ab_50,
    ab_50_severe = ab_50_severe,
    k = k,
    correlated = correlated
  )
  return(parameters)
}


#' @title Combine and verify vaccine parameters
#' @param safir_parameters a list from \code{\link{get_parameters}}
#' @param vaccine_ab_parameters a list from \code{\link{get_vaccine_ab_titre_parameters}}
#' @param vaccine_set a vector giving the number of doses available each day (not each timestep)
#' @param dose_period a vector giving the minimum delay between doses
#' @param strategy_matrix a vaccine strategy matrix from \code{\link[nimue]{strategy_matrix}}
#' or a list of strategy matrices, one for each vaccine dose
#' @param next_dose_priority_matrix a binary matrix giving age groups prioritized for next dose;
#' it should have one fewer row than the number of doses being given, because on the
#' final allocation phase there will be no future dose to prioritize
#' @description Combine parameters for simulation and verify for correctness.
#' @export
make_vaccine_parameters <- function(safir_parameters, vaccine_ab_parameters, vaccine_set, dose_period, strategy_matrix, next_dose_priority_matrix) {

  parameters <- safir_parameters

  vaccine_doses <- length(dose_period)
  stopifnot(vaccine_doses >= 1)
  if (vaccine_doses > 1) {
    stopifnot(all(is.finite(dose_period[2:vaccine_doses])))
    stopifnot(all(dose_period[2:vaccine_doses] >= 0))
  }

  parameters$N_phase <- vaccine_doses
  parameters$dose_period <- dose_period

  stopifnot(length(vaccine_set) == parameters$time_period)
  stopifnot(all(is.finite(vaccine_set)))
  stopifnot(all(vaccine_set >= 0))

  parameters$vaccine_set <- floor(vaccine_set)

  storage.mode(next_dose_priority_matrix) <- "integer"
  stopifnot(nrow(next_dose_priority_matrix) == vaccine_doses - 1)

  # one strategy matrix for all dose phases (input matrix)
  if (inherits(strategy_matrix, 'matrix')) {

    N_age <- ncol(strategy_matrix)
    # N_prioritisation_steps <- nrow(strategy_matrix)

    strategy_matrix_array <- replicate(n = vaccine_doses, expr = strategy_matrix, simplify = FALSE)

  # one strategy matrix for each dose phase (input list)
  } else if(inherits(strategy_matrix, 'list')) {

    stopifnot(length(strategy_matrix) == vaccine_doses)

    N_age <- ncol(strategy_matrix[[1]])
    # N_prioritisation_steps <- nrow(strategy_matrix[[1]])

    strategy_matrix_array <- strategy_matrix

  # bad input
  } else {
    stop("invalid object passed for strategy_matrix")
  }

  stopifnot(ncol(next_dose_priority_matrix) == N_age)
  stopifnot(parameters$N_age == N_age)

  # parameters$N_prioritisation_steps <- N_prioritisation_steps
  parameters$vaccine_coverage_mat <- strategy_matrix_array
  parameters$next_dose_priority <- next_dose_priority_matrix

  parameters <- c(parameters, vaccine_ab_parameters)

  return(parameters)
}
