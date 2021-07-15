# --------------------------------------------------
#   create event functions full vaccine model
#   Sean L. Wu (slwood89@gmail.com)
#   June 2021
# --------------------------------------------------

# write a function

#' @title Append vaccination events
#'
#' @param events a named list of individual::Event
#' @param parameters write me!
#' @export
create_events_vaccination <- function(events, parameters) {

  # pop size
  N <- sum(parameters$population)

  # scheduled future doses
  events$scheduled_dose <- replicate(n = parameters$N_phase,expr = individual::TargetedEvent$new(N),simplify = FALSE)

  return(events)
}


#' @title Schedule some individuals for a vaccination dose
#' @description This is called from the event listeners for each dose, and also
#' aids in better testing of the simulation model
#' @param timestep current time step
#' @param variables a list
#' @param target a \code{\link[individual]{Bitset}}
#' @param dose which dose
#'
#' @export
schedule_dose_vaccine <- function(timestep, variables, target, dose) {

  variables$dose_num$queue_update(value = dose,index = target)
  variables$dose_time[[dose]]$queue_update(values = timestep, index = target)
  if (inherits(target,"Bitset")) {
    variables$ab_titre$queue_update(values = rlnorm(n = target$size()), index = target) # update Ab titre somehow
  } else {
    variables$ab_titre$queue_update(values = rlnorm(n = length(target)), index = target) # update Ab titre somehow
  }


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

  # refactor eventually
  events$exposure$.listeners[[2]] <- NULL
  events$exposure$add_listener(
    create_exposure_scheduler_listener_vaccine(events = events,variables = variables,parameters = parameters,dt = dt)
  )

  # add to each dose
  for (d in seq_along(events$scheduled_dose)) {
    events$scheduled_dose[[d]]$add_listener(
      create_vaccination_dose_listener(variables = variables,parameters = parameters,dose = d)
    )
  }

}
