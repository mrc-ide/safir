# --------------------------------------------------
#   create event functions full vaccine model
#   Sean L. Wu (slwood89@gmail.com)
#   June 2021
# --------------------------------------------------


#' @title Append vaccination events
#'
#' @param events a named list of individual::Event
#' @param parameters write me!
#' @export
create_events_vaccination <- function(events, parameters) {

  # pop size
  N <- sum(parameters$population)

  # scheduled future doses
  events$scheduled_dose <- replicate(n = parameters$N_dose,expr = individual::TargetedEvent$new(N),simplify = FALSE)

  return(events)
}


#' @title Create listener for vaccination dose (multi-dose, no types)
#' @description Updates state when a vaccine dose is given. It does not schedule future events.
#' @param variables a named list of variables
#' @param parameters the parameters
#' @param dose integer dose
create_vaccination_dose_listener <- function(variables, parameters, dose) {
  stopifnot( dose > 0 )
  function (timestep, target) {
    schedule_dose_vaccine(timestep = timestep,variables = variables,target = target,dose = dose)
  }
}


#' @title Attach event listeners for vaccination events (multi-dose, no types)
#' @param variables a named list of variables
#' @param events a named list of events
#' @param parameters the parameters
#' @param type character type of vaccine
#' @export
attach_event_listeners_vaccination <- function(variables, events, parameters, dt) {

  for (d in seq_along(events$scheduled_dose)) {
    events$scheduled_dose[[d]]$add_listener(
      create_vaccination_dose_listener(variables = variables,parameters = parameters,dose = d)
    )
  }

}
