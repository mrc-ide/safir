# --------------------------------------------------
#   exposure event vaccine model with variants of concern (multi-dose, no types)
#   Sean L. Wu (slwood89@gmail.com)
#   September 2021
# --------------------------------------------------


#' @title Modelling the progression to either IMild, ICase, IAsymp with variants of concern (vaccine model, multi-dose, no types)
#' @description Age dependent outcome of exposure
#'
#' @param events a list of events in the model
#' @param variables the available human variables
#' @param parameters model parameters from \code{\link{get_parameters_nimue}}
#' @param dt the time step
#' @param shift passed to \code{\link{make_rerlang}}
#' @export
create_exposure_scheduler_listener_vaccine_voc <- function(events, variables, parameters, dt, shift = 0) {

  stopifnot("voc" %in% names(variables))
  stopifnot("voc_num" %in% names(parameters))

  ICase_delay_voc <- lapply(X = parameters$dur_E, FUN = function(dur_E){make_rerlang(mu = parameters$dur_E, dt = dt, shift = shift)})
  IMild_delay_voc <- lapply(X = parameters$dur_E, FUN = function(dur_E){make_rerlang(mu = parameters$dur_E, dt = dt, shift = shift)})
  IAsymp_delay_voc <- lapply(X = parameters$dur_E, FUN = function(dur_E){make_rerlang(mu = parameters$dur_E, dt = dt, shift = shift)})

  return(
    function(timestep, target) {

      # age of exposed persons
      disc_ages <- variables$discrete_age$get_values(target)

      for (v in seq_len(parameters$voc_num)) {

        # persons exposed with that particular voc
        hosp <- target$copy()
        voc <- variables$voc$get_index_of(parameters$voc_types[v])
        hosp$and(voc)

        if( hosp$size() > 0) {

          # P(hospital | voc, age)
          prob_hosp <- parameters$prob_hosp[v, disc_ages]

          # vaccine efficacy against severe disease
          ab_titre <- variables$ab_titre$get_values(hosp)
          infection_efficacy <- vaccine_efficacy_infection_cpp(ab_titre = ab_titre,parameters = parameters)
          severe_efficacy <- vaccine_efficacy_severe_cpp(ab_titre = ab_titre,ef_infection = infection_efficacy,parameters = parameters)

          # sample those with severe disease
          hosp$sample(prob_hosp * severe_efficacy)

          # those without severe disease
          not_hosp <- target$set_difference(hosp)

          if (hosp$size() > 0) {
            events$severe_infection$schedule(target = hosp, delay = ICase_delay_voc[v](n = hosp$size()))
          }

          # sample asymptomatic and mild disease persons
          if (not_hosp$size() > 0) {
            disc_ages <- variables$discrete_age$get_values(not_hosp)
            prob_asymp <- parameters$prob_asymp[v, disc_ages]

            to_asymp <- not_hosp$copy()
            to_asymp$sample(prob_asymp)

            not_to_asymp <- not_hosp$set_difference(to_asymp)

            if (to_asymp$size() > 0) {
              events$asymp_infection$schedule(target = to_asymp, delay = IAsymp_delay_voc[v](n = to_asymp$size()))
            }

            if (not_to_asymp$size() > 0) {
              events$mild_infection$schedule(target = not_to_asymp, delay = IMild_delay_voc[v](n = not_to_asymp$size()))
            }

          }
        }

      } # end iteration over voc's

    } # end fn
  ) # end return

}

