#' @title Create events
#'
#' @return a named list of individual::Event
#' @noRd
create_events_zzz <- function(parameters, vaccines = FALSE) {

    # pop size
    N <- sum(parameters$population)

    stopifnot(!vaccines)

    list(
        # Human infection events
        exposure = individual::TargetedEvent$new(N), # S->E, scheduled by infection_process_zzz
        mild_infection = individual::TargetedEvent$new(N), # E->IMild, scheduled by create_exposure_update_listener
        asymp_infection = individual::TargetedEvent$new(N),
        severe_infection = individual::TargetedEvent$new(N),
        hospitilisation = individual::TargetedEvent$new(N),
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

#' @title A listener to schedule future events
#' @description a listener to be attached to a \code{\link[individual]{TargetedEvent}}
#' to schedule a future event when it is called.
#' @param event the future event to be schedule
#' @param duration mean duration of waiting time to be scheduled
#' @param func either \code{\link{make_rerlang}} or \code{\link{make_rexp}}
#' @param shift add integer number of time steps to sampled value
#' @param dt size of time step
#' @noRd
create_update_scheduler <- function(event, duration, func, shift, dt) {
  dwell <- func(mu = duration, dt = dt, shift = shift)
  function(timestep, target) {
    event$schedule(target = target, delay = dwell(n = target$size()))
  }
}




#' @title Attach listeners to events
#' @description defines processes for events that can be scheduled in the future
#'
#' @param variables list of variables in the model
#' @param events a list of events in the model
#' @param parameters the model parameters
#' @noRd
attach_event_listeners_zzz <- function(
  variables,
  events,
  parameters
) {

    # Exposure ----------

    # state update
    events$exposure$add_listener(
        create_infection_update_listener(
        variables,
        "E"
        )
    )

    # event scheduling
    events$exposure$add_listener(
        create_exposure_update_listener(
            events,
            variables,
            parameters
        )
    )

    # Mild Infection ----------

    # state update
    events$mild_infection$add_listener(
        create_infection_update_listener(
            variables,
            "IMild"
        )
    )

    # Mild Infection events

    # event scheduling
    events$mild_infection$add_listener(
        create_progression_listener(
            event = events$recovery,
            duration = parameters$dur_IMild,
            func = r_exp
        )
    )



}
