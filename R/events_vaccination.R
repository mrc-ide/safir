# --------------------------------------------------
#   create event functions full vaccine model
#   Sean L. Wu (slwood89@gmail.com)
#   June 2021
# --------------------------------------------------

# write a function

#' @title Append vaccination events
#'
#' @param events a named list of individual::Event
#' @param parameters model parameters
#' @importFrom individual TargetedEvent
#' @export
create_events_vaccination <- function(events, parameters) {

  # pop size
  N <- sum(parameters$population)

  # scheduled future doses
  events$scheduled_dose <- replicate(n = parameters$N_phase,expr = TargetedEvent$new(N),simplify = FALSE)

  return(events)
}


#' @title Schedule some individuals for a vaccination dose
#' @description This is called from the event listeners for each dose, and also
#' aids in better testing of the simulation model
#' @param timestep current time step
#' @param variables a list
#' @param target a \code{\link[individual]{Bitset}}
#' @param dose which dose
#' @param parameters model parameters
#' @importFrom stats rnorm
#' @export
schedule_dose_vaccine <- function(timestep, variables, target, dose, parameters) {

  variables$dose_num$queue_update(value = dose,index = target)
  variables$dose_time[[dose]]$queue_update(values = timestep, index = target)

  if (inherits(target,"Bitset")) {
    n <- target$size()
  } else {
    n <- length(target)
  }

  if (parameters$correlated) {

    if (dose > 1) {
      # correlated doses > 1; use ratio of mean titre
      zdose_prev <- variables$zdose$get_values(index = target)
      zdose_prev <- log10(exp(zdose_prev)) # transform back to log scale
      zdose <- log(10^zdose_prev * (parameters$mu_ab[dose] / parameters$mu_ab[dose-1]))
      variables$zdose$queue_update(values = zdose, index = target)
    } else {
      # initial dose
      zdose <- log(10^rnorm(n = n, mean = log10(parameters$mu_ab[dose]),sd = parameters$std10))
      variables$zdose$queue_update(values = zdose, index = target)
    }

  } else {
    # uncorrelated doses (also use for dose 1 of correlated dose titre)
    zdose <- log(10^rnorm(n = n, mean = log10(parameters$mu_ab[dose]),sd = parameters$std10))
  }

  variables$ab_titre$queue_update(values = zdose, index = target)

}


#' @title Create listener for vaccination dose (multi-dose, no types)
#' @description Updates state when a vaccine dose is given. It does not schedule future events.
#' @param variables a named list of variables
#' @param parameters the parameters
#' @param dose integer dose
create_vaccination_dose_listener <- function(variables, parameters, dose) {
  stopifnot( dose > 0 )
  function (timestep, target) {
    schedule_dose_vaccine(timestep = timestep,variables = variables,target = target,dose = dose, parameters = parameters)
  }
}


#' @title Attach event listeners for vaccination events (multi-dose, no types)
#' @param variables a named list of variables
#' @param events a named list of events
#' @param parameters the parameters
#' @param dt size of time step
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
