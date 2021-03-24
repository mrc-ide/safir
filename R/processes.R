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

#' @title Define event based processes
#' @description defines processes for events that can be scheduled in the future
#'
#' @param human humans
#' @param states a list of states in the model
#' @param variables list of variables in the model
#' @param events a list of events in the model
#' @param parameters the model parameters
#' @noRd
create_event_based_processes <- function(
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
