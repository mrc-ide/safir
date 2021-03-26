#' @title Define event based processes at initialisation
#' @description Defines processes for events for states at initialisation
#'
#' @param parameters model parameters
#' @param events a list of events in the model
#' @param variables list of variables in the model
#' @noRd
create_setup_process <- function(
  parameters,
  events,
  variables
) {
  exposed <- variables$states$get_index_of("E")
  disc_ages <- variables$discrete_age$get_values(exposed)
  prob_hosp <- parameters$prob_hosp[disc_ages]
  hosp <- bernoulli_multi_p(prob_hosp)

  # Get those who have severe infections
  if(sum(hosp) > 0) {
    to_hosp <- individual::filter_bitset(exposed, which(hosp))
    events$severe_infection$schedule(to_hosp,
                                     delay = r_erlang(to_hosp$size(), parameters$dur_E))
  }

  if(sum(!hosp) > 0) {
    # Get individuals not going to hospital
    no_hosp <- which(!hosp)
    not_to_hosp <- individual::filter_bitset(exposed, no_hosp)
    prob_asymp <- parameters$prob_asymp[disc_ages[no_hosp]]
    asymp <- bernoulli_multi_p(prob_asymp)

    # Get those who are asymptomatic
    if (sum(asymp) > 0){
      to_asymp <- individual::filter_bitset(not_to_hosp, which(asymp))
      events$asymp_infection$schedule(to_asymp,
                                      delay = r_erlang(to_asymp$size(),
                                                       parameters$dur_E) )
    }
    # Get those who have mild infections
    if (sum(!asymp) > 0){
      not_to_asymp <- individual::filter_bitset(not_to_hosp, which(!asymp))
      events$mild_infection$schedule(not_to_asymp,
                                     delay = r_erlang(not_to_asymp$size(),
                                                      parameters$dur_E) + 1)
    }

  }
}


#' @title Create list of processes for the simulation
#'
#' @description wires up all the processes to run in the simulation
#'
#' @param human the human handle
#' @param states a list of states in the model
#' @param events a list of events in the model
#' @param variables a list of variables in the model
#' @param parameters a list of parameters in the model
#' @noRd
create_processes <- function(
  events,
  variables,
  parameters,
  renderer
) {

  list(
    infection_process(
      parameters,
      variables,
      events
    ),
    individual::categorical_count_renderer_process(renderer, variables$state, categories = variables$states$get_categories())
  )

}
