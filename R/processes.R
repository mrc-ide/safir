#' @title Define event based processes at initialisation
#' @description Defines processes for events for states at initialisation
#'
#' @param human humans
#' @param states a list of states in the model
#' @param events a list of events in the model
#' @param variables list of variables in the model
#' @noRd
create_setup_process <- function(
  human,
  states,
  events,
  variables
) {
  function(api) {
    parameters <- api$get_parameters()
    exposed <- api$get_state(human, states$E)
    age <- api$get_variable(human, variables$discrete_age, exposed)
    prob_hosp <- parameters$prob_hosp[as.integer(age)]
    hosp <- bernoulli_multi_p(prob_hosp)

    # Get those who have severe infections
    if(sum(hosp) > 0) {
      api$schedule(
        event = events$severe_infection,
        target = exposed[hosp],
        delay = r_erlang(length(exposed[hosp]), parameters$dur_E)
      )
    }

    if(sum(!hosp) > 0) {

      no_hosp <- which(!hosp)
      prob_asymp <-
        parameters$prob_asymp[variables$discrete_age$initial_values[no_hosp]]
      asymp <- bernoulli_multi_p(prob_asymp)

      # Get those who are asymptomatic
      if (sum(asymp) > 0){
        api$schedule(
          event = events$asymp_infection,
          target = exposed[no_hosp][asymp],
          delay = r_erlang(length(exposed[no_hosp][asymp]),
                           parameters$dur_E)
        )
      }
      # Get those who have mild infections
      if (sum(!asymp) > 0){
        api$schedule(
          event = events$mild_infection,
          target = exposed[no_hosp][!asymp],
          delay = r_erlang(length(exposed[no_hosp][!asymp]),
                           parameters$dur_E)
        )

      }

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
  human,
  states,
  events,
  variables,
  parameters
) {

  list(
    infection_process(
      human,
      states,
      variables$discrete_age,
      events$exposure,
      parameters$mix_mat_set
    ),
    individual::state_count_renderer_process(
      human$name,
      unlist(lapply(states, "[[", "name"))
    )
  )

}
