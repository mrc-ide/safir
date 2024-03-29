# --------------------------------------------------
#   exposure event (nimue vaccine model)
#   Sean L. Wu (slwood89@gmail.com)
#   May 2021
#   1. create_exposure_scheduler_listener
# --------------------------------------------------


#' @title Modelling the progression to either IMild, ICase, IAsymp (nimue vaccine model)
#' @description Age dependent outcome of exposure
#'
#' @param events a list of events in the model
#' @param variables the available human variables
#' @param parameters model parameters from \code{\link{get_parameters_nimue}}
#' @param dt the time step
#' @param shift passed to \code{\link{make_rerlang}}
#' @export
create_exposure_scheduler_listener_nimue <- function(events, variables, parameters, dt, shift = 0) {

  stopifnot(length(dim(parameters$prob_hosp)) == 3)
  stopifnot(dim(parameters$prob_hosp)[2] == 17)
  stopifnot(dim(parameters$prob_hosp)[3] == 4)
  stopifnot(all(c("discrete_age", "vaccine_states") %in% names(variables)))

  ICase_delay <- make_rerlang(mu = parameters$dur_E, dt = dt, shift = shift)
  IMild_delay <- make_rerlang(mu = parameters$dur_E, dt = dt, shift = shift)
  IAsymp_delay <- make_rerlang(mu = parameters$dur_E, dt = dt, shift = shift)

  return(
    function(timestep, target) {

      ages <- variables$discrete_age$get_values(target)
      vaxx <- variables$vaccine_states$get_values(target)

      submat <- matrix(data = NA,nrow = target$size(),ncol = 3)
      submat[, 1] <- ceiling(timestep * dt)
      submat[, 2] <- ages
      submat[, 3] <- vaxx

      prob_hosp <- parameters$prob_hosp[submat]
      hosp <- target$copy()

      hosp$sample(prob_hosp)
      not_hosp <- target$set_difference(hosp)

      if (hosp$size() > 0) {
        events$severe_infection$schedule(target = hosp, delay = ICase_delay(n = hosp$size()))
      }

      if (not_hosp$size() > 0) {

        ages <- variables$discrete_age$get_values(not_hosp)

        prob_asymp <- parameters$prob_asymp[ages]

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
