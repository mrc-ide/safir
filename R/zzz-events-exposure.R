#' @title Modelling the progression to either IMild or ICase
#' @description Age dependent outcome of exposure
#'
#' @param events a list of events in the model
#' @param variables the available human variables
#' @param parameters model parameters
#' @noRd
create_exposure_scheduler <- function(events, variables, parameters, dt, vaccines = NULL, shift = 0) {

  ICase_delay <- make_rerlang(mu = parameters$dur_E, dt = dt, shift = shift)
  IMild_delay <- make_rerlang(mu = parameters$dur_E, dt = dt, shift = shift)
  IAsymp_delay <- make_rerlang(mu = parameters$dur_E, dt = dt, shift = shift)

  if (is.null(vaccines)) {

    return(
      function(timestep, to_move) {

        disc_ages <- variables$discrete_age$get_values(to_move)
        prob_hosp <- parameters$prob_hosp[disc_ages]
        hosp <- bernoulli_multi_p(prob_hosp)

        # Severe infections
        if(sum(hosp) > 0) {
          to_hosp <- individual::filter_bitset(to_move, which(hosp))
          events$severe_infection$schedule(target = to_hosp, delay = ICase_delay(n = to_hosp$size()))
        }

        # Non severe infections
        if(sum(!hosp) > 0){
          # Get individuals not going to hospital
          no_hosp <- which(!hosp)
          not_to_hosp <- individual::filter_bitset(to_move, no_hosp)
          prob_asymp <- parameters$prob_asymp[disc_ages[no_hosp]]
          asymp <- bernoulli_multi_p(prob_asymp)

          # Get those who are asymptomatic
          if (sum(asymp) > 0){
            to_asymp <- individual::filter_bitset(not_to_hosp, which(asymp))
            events$asymp_infection$schedule(target = to_asymp, delay = IAsymp_delay(n = to_asymp$size()))

          }
          # Get those who have mild infections
          if (sum(!asymp) > 0){
            not_to_asymp <- individual::filter_bitset(not_to_hosp, which(!asymp))
            events$mild_infection$schedule(target = not_to_asymp, delay = IMild_delay(n = not_to_asymp$size()))

          }

        }
      }
    )

  } else {
    stop("not implemented vaccines yet")
  }

}
