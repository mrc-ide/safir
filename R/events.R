# --------------------------------------------------
#   create event functions
#   May 2021
#   1. create_events
#   2. create_event_scheduler_listener
#   3. create_state_update_listener
#   4. attach_event_listeners
# --------------------------------------------------


#' @title Create events
#'
#' @param parameters model parameters
#' @importFrom individual TargetedEvent
#' @export
create_events <- function(parameters) {

    # pop size
    N <- sum(parameters$population)

    list(
        # Human infection events
        exposure = TargetedEvent$new(N), # S->E, scheduled by infection_process_zzz
        mild_infection = TargetedEvent$new(N), # E->IMild, scheduled by create_exposure_update_listener
        asymp_infection = TargetedEvent$new(N),
        severe_infection = TargetedEvent$new(N),
        hospitilisation = TargetedEvent$new(N),
        imv_get_live = TargetedEvent$new(N),
        imv_get_die = TargetedEvent$new(N),
        iox_get_live = TargetedEvent$new(N),
        iox_get_die = TargetedEvent$new(N),
        imv_not_get_live = TargetedEvent$new(N),
        imv_not_get_die = TargetedEvent$new(N),
        iox_not_get_live = TargetedEvent$new(N),
        iox_not_get_die = TargetedEvent$new(N),
        stepdown = TargetedEvent$new(N),
        recovery = TargetedEvent$new(N),
        immunity_loss = TargetedEvent$new(N),
        death = TargetedEvent$new(N)
    )
}


#' @title A listener to schedule future events
#' @description a listener to be attached to a \code{\link[individual]{TargetedEvent}}
#' to schedule a future event when that event fires.
#' @param event the future event to be schedule
#' @param duration mean duration of waiting time to be scheduled
#' @param func either \code{\link{make_rerlang}} or \code{\link{make_rexp}}
#' @param shift add integer number of time steps to sampled value
#' @param dt size of time step
#' @export
create_event_scheduler_listener <- function(event, duration, func, shift, dt) {
  dwell <- func(mu = duration, dt = dt, shift = shift)
  function(timestep, target) {
    event$schedule(target = target, delay = dwell(n = target$size()))
  }
}


#' @title A listener to update state
#' @description a listener to be attached to a \code{\link[individual]{TargetedEvent}}
#' to update state when that event fires.
#' @param states a \code{\link[individual]{CategoricalVariable}} object
#' @param destination the destination state
#' @export
create_state_update_listener <- function(states, destination) {
    function(timestep, target) {
        states$queue_update(value = destination, index = target)
    }
}


#' @title Attach listeners to events
#' @description defines processes for events that can be scheduled in the future
#'
#' @param variables list of variables in the model
#' @param events a list of events in the model
#' @param parameters the model parameters
#' @param dt size of time step
#' @param shift schedule future events after minimum number of time step delay
#' @param shift_exposure schedule exposure event after minimum number of time step delay
#' @export
attach_event_listeners <- function(
  variables,
  events,
  parameters,
  dt,
  shift = 1L,
  shift_exposure = 1L
) {

    # Exposure ----------

    events$exposure$add_listener(
        create_state_update_listener(
            variables$states,
            "E"
        )
    )

    events$exposure$add_listener(
        create_exposure_scheduler_listener(
            events,
            variables,
            parameters,
            dt = dt,
            shift = shift_exposure
        )
    )

    # IMild ----------

    events$mild_infection$add_listener(
        create_state_update_listener(
            variables$states,
            "IMild"
        )
    )

    events$mild_infection$add_listener(
        create_event_scheduler_listener(
            event = events$recovery,
            duration = parameters$dur_IMild,
            func = make_rexp,
            shift = shift,
            dt = dt
        )
    )

    # IAsymp ----------

    events$asymp_infection$add_listener(
        create_state_update_listener(
            variables$states,
            "IAsymp"
        )
    )

    events$asymp_infection$add_listener(
        create_event_scheduler_listener(
            event = events$recovery,
            duration = parameters$dur_IAsymp,
            func = make_rexp,
            shift = shift,
            dt = dt
        )
    )

    # ICase ----------

    events$severe_infection$add_listener(
        create_state_update_listener(
            variables$states,
            "ICase"
        )
    )

    events$severe_infection$add_listener(
        create_event_scheduler_listener(
            event = events$hospitilisation,
            duration = parameters$dur_ICase,
            func = make_rerlang,
            shift = shift,
            dt = dt
        )
    )

    # Hospitalisation (no state update, queues other events) ----------

    events$hospitilisation$add_listener(
        create_hospital_scheduler_listener(
            parameters = parameters,
            variables = variables,
            events = events
        )
    )

    # IMV (hospitalised, mechanical ventilation) ----------

    events$imv_get_live$add_listener(
        create_state_update_listener(
            variables$states,
            "IMVGetLive"
        )
    )

    events$imv_get_live$add_listener(
        create_event_scheduler_listener(
            event = events$stepdown,
            duration = parameters$dur_get_mv_survive,
            func = make_rerlang,
            shift = shift,
            dt = dt
        )
    )

    events$imv_get_die$add_listener(
        create_state_update_listener(
            variables$states,
            "IMVGetDie"
        )
    )

    events$imv_get_die$add_listener(
        create_event_scheduler_listener(
            event = events$death,
            duration = parameters$dur_get_mv_die,
            func = make_rerlang,
            shift = shift,
            dt = dt
        )
    )

    events$imv_not_get_live$add_listener(
        create_state_update_listener(
            variables$states,
            "IMVNotGetLive"
        )
    )

    events$imv_not_get_live$add_listener(
        create_event_scheduler_listener(
            event = events$recovery,
            duration = parameters$dur_not_get_mv_survive,
            func = make_rerlang,
            shift = shift,
            dt = dt
        )
    )

    events$imv_not_get_die$add_listener(
        create_state_update_listener(
            variables$states,
            "IMVNotGetDie"
        )
    )

    events$imv_not_get_die$add_listener(
        create_event_scheduler_listener(
            event = events$death,
            duration = parameters$dur_not_get_mv_die,
            func = make_rerlang,
            shift = shift,
            dt = dt
        )
    )

    # IOx (hospitalised, oxygen) ----------

    events$iox_get_live$add_listener(
        create_state_update_listener(
            variables$states,
            "IOxGetLive"
        )
    )

    events$iox_get_live$add_listener(
        create_event_scheduler_listener(
            event = events$recovery,
            duration = parameters$dur_get_ox_survive,
            func = make_rerlang,
            shift = shift,
            dt = dt
        )
    )

    events$iox_get_die$add_listener(
        create_state_update_listener(
            variables$states,
            "IOxGetDie"
        )
    )

    events$iox_get_die$add_listener(
        create_event_scheduler_listener(
            event = events$death,
            duration = parameters$dur_get_ox_die,
            func = make_rerlang,
            shift = shift,
            dt = dt
        )
    )

    events$iox_not_get_live$add_listener(
        create_state_update_listener(
            variables$states,
            "IOxNotGetLive"
        )
    )

    events$iox_not_get_live$add_listener(
        create_event_scheduler_listener(
            event = events$recovery,
            duration = parameters$dur_not_get_ox_survive,
            func = make_rerlang,
            shift = shift,
            dt = dt
        )
    )

    events$iox_not_get_die$add_listener(
        create_state_update_listener(
            variables$states,
            "IOxNotGetDie"
        )
    )

    events$iox_not_get_die$add_listener(
        create_event_scheduler_listener(
            event = events$death,
            duration = parameters$dur_not_get_ox_die,
            func = make_rerlang,
            shift = shift,
            dt = dt
        )
    )


    # Recovery events
    events$recovery$add_listener(
        create_state_update_listener(
            variables$states,
            "R"
        )
    )

    if (is.finite(parameters$dur_R)) {
        events$recovery$add_listener(
            create_event_scheduler_listener(
                event = events$immunity_loss,
                duration = parameters$dur_R,
                func = make_rerlang,
                shift = shift,
                dt = dt
            )
        )
    }

    # Stepdown events
    events$stepdown$add_listener(
        create_state_update_listener(
            variables$states,
            "IRec"
        )
    )

    events$stepdown$add_listener(
        create_event_scheduler_listener(
            event = events$recovery,
            duration = parameters$dur_rec,
            func = make_rerlang,
            shift = shift,
            dt = dt
        )
    )

    # Death events
    events$death$add_listener(
        create_state_update_listener(
            variables$states,
            "D"
        )
    )

    # Loss of immunity
    events$immunity_loss$add_listener(
        create_state_update_listener(
            variables$states,
            "S"
        )
    )



}
