#' @ Create events
#'
#' @return a named list of individual::Event
#' @noRd
create_events <- function() {
  list(
    # Human infection events
    exposure = individual::Event$new('exposure'),
    mild_infection = individual::Event$new('mild_infection'),
    asymp_infection = individual::Event$new('asymp_infection'),
    severe_infection = individual::Event$new('severe_infection'), # requiring hospital eventually
    hospitilisation = individual::Event$new('hospitilisation'), # either ICU or MV
    imv_get_live = individual::Event$new('imv_get_live'),
    imv_get_die = individual::Event$new('imv_get_die'),
    iox_get_live = individual::Event$new('iox_get_live'),
    iox_get_die = individual::Event$new('iox_get_die'),
    imv_not_get_live = individual::Event$new('imv_not_get_live'),
    imv_not_get_die = individual::Event$new('imv_not_get_die'),
    iox_not_get_live = individual::Event$new('iox_not_get_live'),
    iox_not_get_die = individual::Event$new('iox_not_get_die'),
    stepdown = individual::Event$new('stepdown'),
    recovery = individual::Event$new('recovery'),
    immunity_loss = individual::Event$new('immunity_loss'),
    death = individual::Event$new('death')
  )
}

#' @title Update the state of an individual as infection events occur
#' @description Moves individuals towards the later stages of disease
#'
#' @param human the handle for the human individuals
#' @param to_state the destination disease state
#' @noRd
create_infection_update_listener <- function(
  human,
  to_state) {
  function(api, to_move) {
    api$queue_state_update(human, to_state, to_move)
  }
}

#' @title Schedule progression of human disease at the start of the simulation
#' @description Schedules infection events using Erlang
#'
#' @param event the event to schedule
#' @param human the human handle
#' @param from_state the state this event applies to
#' @param duration the average time spent in this state
#' @noRd
initialise_progression <- function(event, human, from_state, duration) {
  function(api, target) {
    target <- api$get_state(human, from_state)
    api$schedule(event, target, r_erlang(length(target), duration))
  }
}

#' @title Modelling the progression of the human disease
#' @description schedules follow up infection events
#'
#' @param event the event to schedule
#' @param duration the average time spent in this state
#' @param shift number of days to increase scheduled event by
#' @param func function to use for drawing progression.
#' Default = [r_erlang]
#' @noRd
create_progression_listener <- function(event, duration, shift = 0, func = r_erlang) {
  function(api, target) {
    api$schedule(event, target, func(length(target), duration) + shift)
  }
}

#' @title Modelling the progression to either IMild or ICase
#' @description Age dependent outcome of exposure
#'
#' @param human the handle for the human individuals
#' @param states the available human states
#' @param events a list of events in the model
#' @param variables the available human variables
#' @param parameters model parameters
#' @noRd
create_exposure_update_listener <- function(
  human,
  states,
  events,
  variables,
  parameters) {
  function(api, to_move) {
    disc_ages <- api$get_variable(human, variables$discrete_age, to_move)
    prob_hosp <- parameters$prob_hosp[as.integer(disc_ages)]
    hosp <- bernoulli_multi_p(prob_hosp)

    # Severe infections
    if(sum(hosp) > 0) {
      api$schedule(
        event = events$severe_infection,
        target = to_move[hosp],
        delay = r_erlang(length(to_move[hosp]), parameters$dur_E) + 1
      )
    }

    # Non severe infections
    if(sum(!hosp) > 0){
      # Get individuals not going to hospital
      no_hosp <- which(!hosp)
      prob_asymp <- parameters$prob_asymp[as.integer(disc_ages[no_hosp])]
      asymp <- bernoulli_multi_p(prob_asymp)

      # Get those who are asymptomatic
      if (sum(asymp) > 0){
        api$schedule(
          event = events$asymp_infection,
          target = to_move[no_hosp][asymp],
          delay = r_erlang(length(to_move[no_hosp][asymp]),
                           parameters$dur_E) + 1
        )
      }
      # Get those who have mild infections
      if (sum(!asymp) > 0){
        api$schedule(
          event = events$mild_infection,
          target = to_move[no_hosp][!asymp],
          delay = r_erlang(length(to_move[no_hosp][!asymp]),
                           parameters$dur_E) + 1
        )
      }

    }
  }
}

