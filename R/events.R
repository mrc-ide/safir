#' @ Create events
#'
#' @return a named list of individual::Event
#' @noRd
create_events <- function(parameters) {

  # pop size
  N <- sum(parameters$population)

  list(
    # Human infection events
    exposure = individual::TargetedEvent$new(N),
    mild_infection = individual::TargetedEvent$new(N),
    asymp_infection = individual::TargetedEvent$new(N),
    severe_infection = individual::TargetedEvent$new(N), # requiring hospital eventually
    hospitilisation = individual::TargetedEvent$new(N), # either ICU or MV
    imv_get_live = individual::TargetedEvent$new(N),
    imv_get_die = individual::TargetedEvent$new(N),
    iox_get_live = individual::TargetedEvent$new(N),
    iox_get_die = individual::TargetedEvent$new(N),
    imv_not_get_live = individual::TargetedEvent$new(N),
    imv_not_get_die = individual::TargetedEvent$new(N),
    iox_not_get_live = individual::TargetedEvent$new(N),
    iox_not_get_die = individual::TargetedEvent$new(N),
    stepdown = individual::TargetedEvent$new(N),
    recovery = individual::TargetedEvent$new(N),
    immunity_loss = individual::TargetedEvent$new(N),
    death = individual::TargetedEvent$new(N)
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


#' @title Attach listeners to events
#' @description defines processes for events that can be scheduled in the future
#'
#' @param human humans
#' @param states a list of states in the model
#' @param variables list of variables in the model
#' @param events a list of events in the model
#' @param parameters the model parameters
#' @noRd
attach_event_listeners <- function(
  human,
  states,
  variables,
  events,
  parameters
) {

  # STATE UPDATES
  # These events cause the infection state to change at the end of the timestep
  # ---------------------
  # Exposure events
  events$exposure$add_listener(
    create_infection_update_listener(
      human,
      states$E
    )
  )

  # IMild events
  events$mild_infection$add_listener(
    create_infection_update_listener(
      human,
      states$IMild
    )
  )

  # IAsymp events
  events$asymp_infection$add_listener(
    create_infection_update_listener(
      human,
      states$IAsymp
    )
  )

  # ICase events
  events$severe_infection$add_listener(
    create_infection_update_listener(
      human,
      states$ICase
    )
  )

  # IMV events
  events$imv_get_live$add_listener(
    create_infection_update_listener(
      human,
      states$IMVGetLive
    )
  )

  events$imv_get_die$add_listener(
    create_infection_update_listener(
      human,
      states$IMVGetDie
    )
  )

  events$imv_not_get_live$add_listener(
    create_infection_update_listener(
      human,
      states$IMVNotGetLive
    )
  )

  events$imv_not_get_die$add_listener(
    create_infection_update_listener(
      human,
      states$IMVNotGetDie
    )
  )

  # IOx events
  events$iox_get_live$add_listener(
    create_infection_update_listener(
      human,
      states$IOxGetLive
    )
  )

  events$iox_get_die$add_listener(
    create_infection_update_listener(
      human,
      states$IOxGetDie
    )
  )

  events$iox_not_get_live$add_listener(
    create_infection_update_listener(
      human,
      states$IOxNotGetLive
    )
  )

  events$iox_not_get_die$add_listener(
    create_infection_update_listener(
      human,
      states$IOxNotGetDie
    )
  )

  # Recovery events
  events$recovery$add_listener(
    create_infection_update_listener(
      human,
      states$R
    )
  )

  # Stepdown events
  events$stepdown$add_listener(
    create_infection_update_listener(
      human,
      states$IRec
    )
  )

  # Death events
  events$death$add_listener(
    create_infection_update_listener(
      human,
      states$D
    )
  )

  # Loss of immunity
  events$immunity_loss$add_listener(
    create_infection_update_listener(
      human,
      states$S
    )
  )

  # STATE SCHEDULES
  # These events trigger the scheduling for infection state changes
  # ----------------------------

  # Exposure events
  events$exposure$add_listener(
    create_exposure_update_listener(
      human,
      states,
      events,
      variables,
      parameters
    )
  )

  # Mild Infection events
  events$mild_infection$add_listener(
    create_progression_listener(
      event = events$recovery,
      duration = parameters$dur_IMild,
      func = r_exp
    )
  )

  # Asymptomatic Infection events
  events$asymp_infection$add_listener(
    create_progression_listener(
      event = events$recovery,
      duration = parameters$dur_IAsymp,
      func = r_exp
    )
  )

  # Severe Infection events
  events$severe_infection$add_listener(
    create_progression_listener(
      event = events$hospitilisation,
      duration = parameters$dur_ICase
    )
  )

  # Hospitalisation
  events$hospitilisation$add_listener(
    hospitilisation_flow_process(
      variables$discrete_age,
      human,
      states,
      events
    )
  )

  # MV outcomes
  events$imv_get_live$add_listener(
    create_progression_listener(
      event = events$stepdown,
      duration = parameters$dur_get_mv_survive,
      shift = 1
    )
  )

  events$imv_get_die$add_listener(
    create_progression_listener(
      event = events$death,
      duration = parameters$dur_get_mv_die,
      shift = 1
    )
  )

  events$imv_not_get_live$add_listener(
    create_progression_listener(
      event = events$recovery,
      duration = parameters$dur_not_get_mv_survive,
      shift = 1
    )
  )

  events$imv_not_get_die$add_listener(
    create_progression_listener(
      event = events$death,
      duration = parameters$dur_not_get_mv_die,
      shift = 1
    )
  )

  # Ox outcomes
  events$iox_get_live$add_listener(
    create_progression_listener(
      event = events$recovery,
      duration = parameters$dur_get_ox_survive,
      shift = 1
    )
  )

  events$iox_get_die$add_listener(
    create_progression_listener(
      event = events$death,
      duration = parameters$dur_get_ox_die,
      shift = 1
    )
  )

  events$iox_not_get_live$add_listener(
    create_progression_listener(
      event = events$recovery,
      duration = parameters$dur_not_get_ox_survive,
      shift = 1
    )
  )

  events$iox_not_get_die$add_listener(
    create_progression_listener(
      event = events$death,
      duration = parameters$dur_not_get_ox_die,
      shift = 1
    )
  )

  # Stepdown
  events$stepdown$add_listener(
    create_progression_listener(
      event = events$recovery,
      duration = parameters$dur_rec
    )
  )

  # Loss of Immunity
  if (is.finite(parameters$dur_R)) {
    events$recovery$add_listener(
      create_progression_listener(
        event = events$immunity_loss,
        duration = parameters$dur_R
      )
    )
  }

}
