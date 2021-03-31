#' @title Create events
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
#' @param variables the handle for the variables
#' @param to_state the destination disease state
#' @noRd
create_infection_update_listener <- function(
  variables,
  to_state) {
  function(timestep, to_move) {
    variables$states$queue_update(to_state, to_move)
  }
}

#' @title Schedule progression of human disease at the start of the simulation
#' @description Schedules infection events using Erlang
#'
#' @param variables the handle for the variables
#' @param event the event to schedule
#' @param from_state the state this event applies to
#' @param duration the average time spent in this state
#' @noRd
initialise_progression <- function(variables, event, from_state, duration) {
  function(timestep) {
    target <- variables$state$get_index_of(from_state)
    event$schedule(target, r_erlang(target$size(), duration))
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
  function(timestep, target) {
    event$schedule(target, func(target$size(), duration) + shift)
  }
}

#' @title Modelling the progression to either IMild or ICase
#' @description Age dependent outcome of exposure
#'
#' @param events a list of events in the model
#' @param variables the available human variables
#' @param parameters model parameters
#' @noRd
create_exposure_update_listener <- function(
  events,
  variables,
  parameters) {
  function(timestep, to_move) {
    disc_ages <- variables$discrete_age$get_values(to_move)
    prob_hosp <- parameters$prob_hosp[disc_ages]
    hosp <- bernoulli_multi_p(prob_hosp)

    # Severe infections
    if(sum(hosp) > 0) {
      to_hosp <- individual::filter_bitset(to_move, which(hosp))
      events$severe_infection$schedule(to_hosp,
                                       delay = r_erlang(to_hosp$size(), parameters$dur_E) + 1)
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
        events$asymp_infection$schedule(to_asymp,
                                        delay = r_erlang(to_asymp$size(),
                                                         parameters$dur_E) + 1)

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
}


#' @title Attach listeners to events
#' @description defines processes for events that can be scheduled in the future
#'
#' @param variables list of variables in the model
#' @param events a list of events in the model
#' @param parameters the model parameters
#' @noRd
attach_event_listeners <- function(
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
      variables,
      "E"
    )
  )

  # IMild events
  events$mild_infection$add_listener(
    create_infection_update_listener(
      variables,
      "IMild"
    )
  )

  # IAsymp events
  events$asymp_infection$add_listener(
    create_infection_update_listener(
      variables,
      "IAsymp"
    )
  )

  # ICase events
  events$severe_infection$add_listener(
    create_infection_update_listener(
      variables,
      "ICase"
    )
  )

  # IMV events
  events$imv_get_live$add_listener(
    create_infection_update_listener(
      variables,
      "IMVGetLive"
    )
  )

  events$imv_get_die$add_listener(
    create_infection_update_listener(
      variables,
      "IMVGetDie"
    )
  )

  events$imv_not_get_live$add_listener(
    create_infection_update_listener(
      variables,
      "IMVNotGetLive"
    )
  )

  events$imv_not_get_die$add_listener(
    create_infection_update_listener(
      variables,
      "IMVNotGetDie"
    )
  )

  # IOx events
  events$iox_get_live$add_listener(
    create_infection_update_listener(
      variables,
      "IOxGetLive"
    )
  )

  events$iox_get_die$add_listener(
    create_infection_update_listener(
      variables,
      "IOxGetDie"
    )
  )

  events$iox_not_get_live$add_listener(
    create_infection_update_listener(
      variables,
      "IOxNotGetLive"
    )
  )

  events$iox_not_get_die$add_listener(
    create_infection_update_listener(
      variables,
      "IOxNotGetDie"
    )
  )

  # Recovery events
  events$recovery$add_listener(
    create_infection_update_listener(
      variables,
      "R"
    )
  )

  # Stepdown events
  events$stepdown$add_listener(
    create_infection_update_listener(
      variables,
      "IRec"
    )
  )

  # Death events
  events$death$add_listener(
    create_infection_update_listener(
      variables,
      "D"
    )
  )

  # Loss of immunity
  events$immunity_loss$add_listener(
    create_infection_update_listener(
      variables,
      "S"
    )
  )

  # STATE SCHEDULES
  # These events trigger the scheduling for infection state changes
  # ----------------------------

  # Exposure events
  events$exposure$add_listener(
    create_exposure_update_listener(
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
      parameters,
      variables,
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
