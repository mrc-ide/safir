# --------------------------------------------------
#   create event functions
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
#' @export
create_events_nimue <- function(events, parameters) {

  # pop size
  N <- sum(parameters$population)

  events$v0_to_v1v2 = individual::TargetedEvent$new(N) # scheduled when vaccination occurs
  events$v1v2_to_v3v4 = individual::TargetedEvent$new(N) # scheduled when entering v1v2
  events$v3v4_to_v5 = individual::TargetedEvent$new(N) # scheduled when entering v3v4

  return(events)
}


#' @title Create listener for initial vaccination (nimue vaccine model)
#'
#' @param variables a named list of variables
#' @param events a named list of individual::Event
#' @param parameters the parameters
#' @param func function to draw waiting times from
#' @param shift number of time steps to shift scheduled event by
#' @param dt size of time step
create_v0_to_v1v2_listener_nimue <- function(variables, events, parameters, func, shift, dt) {

  dwell <- func(mu = parameters$gamma_vaccine_delay, dt = dt, shift = shift)

  stopifnot(all(nimue_events_names %in% names(events)))
  stopifnot(all(nimue_states_names %in% names(variables)))

  function (timestep, target) {
    # newly vaccinated persons: update state
    variables$vaccinated$or(target)
    variables$vaccine_states$queue_update(values = 2, index = target)
    # schedule future events
    event$v1v2_to_v3v4$schedule(target = target, delay = dwell(n = target$size()))
  }
}


#' @title Create listener for start of vaccine protection (nimue vaccine model)
#'
#' @param variables a named list of variables
#' @param events a named list of individual::Event
#' @param parameters the parameters
#' @param func function to draw waiting times from
#' @param shift number of time steps to shift scheduled event by
#' @param dt size of time step
create_v1v2_to_v3v4_listener_nimue <- function(variables, events, parameters, func, shift, dt) {

  dwell <- func(mu = parameters$gamma_V, dt = dt, shift = shift)

  stopifnot(all(nimue_events_names %in% names(events)))
  stopifnot(all(nimue_states_names %in% names(variables)))

  function (timestep, target) {
    # newly protected persons: update state
    variables$vaccine_states$queue_update(values = 3, index = target)
    # schedule future events
    event$v3v4_to_v5$schedule(target = target, delay = dwell(n = target$size()))
  }
}


#' @title Create listener for decay of vaccine protection (nimue vaccine model)
#'
#' @param variables a named list of variables
#' @param events a named list of individual::Event
#' @param parameters the parameters
#' @param func function to draw waiting times from
#' @param shift number of time steps to shift scheduled event by
#' @param dt size of time step
create_v3v4_to_v5_listener_nimue <- function(variables, events, parameters, func, shift, dt) {

  stopifnot(all(nimue_events_names %in% names(events)))
  stopifnot(all(nimue_states_names %in% names(variables)))

  function (timestep, target) {
    # newly decay of protection persons: update state
    variables$vaccine_states$queue_update(values = 4, index = target)
  }
}



attach_listeners_nimue <- function() {}

