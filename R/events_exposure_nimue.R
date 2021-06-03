# --------------------------------------------------
#   exposure event
#   May 2021
#   1. create_exposure_scheduler_listener
# --------------------------------------------------


#' @title Modelling the progression to either IMild or ICase
#' @description Age dependent outcome of exposure
#'
#' @param events a list of events in the model
#' @param variables the available human variables
#' @param parameters model parameters
#' @param dt the time step
#' @param shift passed to \code{\link{make_rerlang}}
#' @noRd
create_exposure_scheduler_listener_nimue <- function(events, variables, parameters, dt, shift = 0) {

  ICase_delay <- make_rerlang(mu = parameters$dur_E, dt = dt, shift = shift)
  IMild_delay <- make_rerlang(mu = parameters$dur_E, dt = dt, shift = shift)
  IAsymp_delay <- make_rerlang(mu = parameters$dur_E, dt = dt, shift = shift)

  return(
    function(timestep, to_move) {

      disc_ages <- variables$discrete_age$get_values(to_move)
      prob_hosp <- parameters$prob_hosp[disc_ages]
      hosp <- to_move$copy()

      hosp$sample(prob_hosp)
      not_hosp <- to_move$set_difference(hosp)

      if (hosp$size() > 0) {
        events$severe_infection$schedule(target = hosp, delay = ICase_delay(n = hosp$size()))
      }

      if (not_hosp$size() > 0) {
        disc_ages <- variables$discrete_age$get_values(not_hosp)
        prob_asymp <- parameters$prob_asymp[disc_ages]

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
