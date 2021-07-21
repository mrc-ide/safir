# --------------------------------------------------
#   create event functions (nimue vaccine model)
#   Sean L. Wu (slwood89@gmail.com)
#   May 2021
#   1. create_events
#   2. create_event_scheduler_listener
#   3. create_state_update_listener
#   4. attach_event_listeners
# --------------------------------------------------

nimue_events_names <- c("v0_to_v1v2", "v1v2_to_v3v4", "v3v4_to_v5")
nimue_states_names <- c("vaccine_states", "eligible", "vaccinated", "empty")


#' @title Append vaccination events (nimue vaccine model)
#'
#' @param events a named list of individual::Event
#' @param parameters write me!
#' @importFrom individual TargetedEvent
#' @export
create_events_nimue <- function(events, parameters) {

  # pop size
  N <- sum(parameters$population)

  events$v0_to_v1v2 = TargetedEvent$new(N) # scheduled when vaccination occurs
  events$v1v2_to_v3v4 = TargetedEvent$new(N) # scheduled when entering v1v2
  events$v3v4_to_v5 = TargetedEvent$new(N) # scheduled when entering v3v4

  return(events)
}


#' @title Create listener for initial vaccination (nimue vaccine model)
#' @description This is called when entering v1v2 state
#' @param variables a named list of variables
#' @param events a named list of individual::Event
#' @param parameters the parameters
#' @param func function to draw waiting times from
#' @param shift number of time steps to shift scheduled event by
#' @param dt size of time step
create_v0_to_v1v2_listener_nimue <- function(variables, events, parameters, func, shift, dt) {

  stopifnot(all(nimue_events_names %in% names(events)))
  stopifnot(all(nimue_states_names %in% names(variables)))

  dwell <- func(mu = parameters$dur_vaccine_delay, dt = dt, shift = shift)

  function (timestep, target) {
    # newly vaccinated persons: update state
    variables$vaccinated$or(target)
    variables$vaccine_states$queue_update(values = 2, index = target)
    # schedule future events
    events$v1v2_to_v3v4$schedule(target = target, delay = dwell(n = target$size()))
  }
}


#' @title Create listener for start of vaccine protection (nimue vaccine model)
#' @description This is called when entering v3v4 state
#' @param variables a named list of variables
#' @param events a named list of individual::Event
#' @param parameters the parameters
#' @param func function to draw waiting times from
#' @param shift number of time steps to shift scheduled event by
#' @param dt size of time step
create_v1v2_to_v3v4_listener_nimue <- function(variables, events, parameters, func, shift, dt) {

  stopifnot(all(nimue_events_names %in% names(events)))
  stopifnot(all(nimue_states_names %in% names(variables)))

  dwell <- func(mu = parameters$dur_V, dt = dt, shift = shift)

  function (timestep, target) {
    # newly protected persons: update state
    variables$vaccine_states$queue_update(values = 3, index = target)
    # schedule future events
    events$v3v4_to_v5$schedule(target = target, delay = dwell(n = target$size()))
  }
}


#' @title Create listener for decay of vaccine protection (nimue vaccine model)
#' @description This is called when entering v5 state
#' @param variables a named list of variables
#' @param events a named list of individual::Event
#' @param parameters the parameters
create_v3v4_to_v5_listener_nimue <- function(variables, events, parameters) {

  stopifnot(all(nimue_events_names %in% names(events)))
  stopifnot(all(nimue_states_names %in% names(variables)))

  function (timestep, target) {
    # newly decay of protection persons: update state
    variables$vaccine_states$queue_update(values = 4, index = target)
  }
}


#' @title Attach listeners to events (nimue vaccine model)
#' @description defines processes for events that can be scheduled in the future
#'
#' @param variables list of variables in the model
#' @param events a list of events in the model
#' @param parameters the model parameters
#' @param dt size of time step
#' @param shift schedule future events after minimum number of time step delay
#' @export
attach_event_listeners_nimue <- function(
  variables,
  events,
  parameters,
  dt,
  shift = 1L
) {

  # vaccination (v0 -> v1,v2) ------------------------------
  events$v0_to_v1v2$add_listener(
    create_v0_to_v1v2_listener_nimue(
      variables = variables,
      events = events,
      parameters = parameters,
      func = make_rerlang,
      shift = shift,
      dt = dt
    )
  )

  # start of vaccine protection (v1,v2 -> v3,v4) ------------------------------
  events$v1v2_to_v3v4$add_listener(
    create_v1v2_to_v3v4_listener_nimue(
      variables = variables,
      events = events,
      parameters = parameters,
      func = make_rerlang,
      shift = shift,
      dt = dt
    )
  )

  # decay of vaccine protection (v3,v4 -> v5) ------------------------------
  events$v3v4_to_v5$add_listener(
    create_v3v4_to_v5_listener_nimue(
      variables = variables,
      events = events,
      parameters = parameters
    )
  )

  # modify exposure event
  events$exposure$.listeners[[2]] <- NULL
  events$exposure$add_listener(
    create_exposure_scheduler_listener_nimue(
      events = events,
      variables = variables,
      parameters = parameters,
      dt = dt,
      shift = 0L
    )
  )

}

