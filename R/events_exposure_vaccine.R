# --------------------------------------------------
#   exposure event vaccine model (multi-dose, no types)
#   Sean L. Wu (slwood89@gmail.com)
#   July 2021
# --------------------------------------------------


#' @title Modelling the progression to either IMild, ICase, IAsymp (vaccine model, multi-dose, no types)
#' @description Age dependent outcome of exposure
#'
#' @param events a list of events in the model
#' @param variables the available human variables
#' @param parameters model parameters from \code{\link{get_parameters_nimue}}
#' @param dt the time step
#' @param shift passed to \code{\link{make_rerlang}}
#' @export
create_exposure_scheduler_listener_vaccine <- function(events, variables, parameters, dt, shift = 0) {

  ICase_delay <- make_rerlang(mu = parameters$dur_E, dt = dt, shift = shift)
  IMild_delay <- make_rerlang(mu = parameters$dur_E, dt = dt, shift = shift)
  IAsymp_delay <- make_rerlang(mu = parameters$dur_E, dt = dt, shift = shift)

  calculate_nat <- make_calculate_nat(variables = variables)

  check_probabilities(prob = parameters$prob_hosp, parameters = parameters)
  check_probabilities(prob = parameters$prob_asymp, parameters = parameters)

  return(
    function(timestep, target) {

      # current day of simulation
      day <- ceiling(timestep * dt)

      # probabilities of hospitalization by age group
      disc_ages <- variables$discrete_age$get_values(target)
      # prob_hosp <- parameters$prob_hosp[disc_ages]
      prob_hosp <- get_probabilties(prob = parameters$prob_hosp, ages = disc_ages, day = day)

      hosp <- target$copy()

      # vaccine efficacy against severe disease
      nat <- calculate_nat(variables = variables, index = hosp)
      infection_efficacy <- vaccine_efficacy_infection_cpp(ab_titre = nat, parameters = parameters, day = day - 1L) # 0-based index
      severe_efficacy <- vaccine_efficacy_severe_cpp(ab_titre = nat, ef_infection = infection_efficacy,parameters = parameters, day = day - 1L) # 0-based index

      # sample those with severe disease
      hosp$sample(prob_hosp * severe_efficacy)

      # those without severe disease
      not_hosp <- target$set_difference(hosp)

      if (hosp$size() > 0) {
        events$severe_infection$schedule(target = hosp, delay = ICase_delay(n = hosp$size()))
      }

      # sample asymptomatic and mild disease persons
      if (not_hosp$size() > 0) {
        disc_ages <- variables$discrete_age$get_values(not_hosp)
        # prob_asymp <- parameters$prob_asymp[disc_ages]
        prob_asymp <- get_probabilties(prob = parameters$prob_asymp, ages = disc_ages, day = day)

        to_asymp <- not_hosp$copy()
        to_asymp$sample(prob_asymp)

        not_to_asymp <- not_hosp$set_difference(to_asymp)

        if (to_asymp$size() > 0) {
          events$asymp_infection$schedule(target = to_asymp, delay = IAsymp_delay(n = to_asymp$size()))
        }

        if (not_to_asymp$size() > 0) {
          events$mild_infection$schedule(target = not_to_asymp, delay = IMild_delay(n = not_to_asymp$size()))
        }

      }

    }
  )

}
