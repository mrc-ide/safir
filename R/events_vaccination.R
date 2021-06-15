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
#' @param type character type of vaccine
#' @param dose integer dose
create_vaccination_dose_listener <- function(variables, parameters, type, dose) {
  stopifnot(type %in% parameters$vaxx_types)
  stopifnot( parameters$N_dose[which(parameters$vaxx_types == type)] >= dose )
  stopifnot( dose > 0 )
  function (timestep, target) {
    variables$dose_num$queue_update(values = dose, index = target)
    variables$dose_time[[dose]]$queue_update(values = timestep, index = target)
  }
}


# types

# create_events_vaccination <- function(events, parameters) {
#
#   # pop size
#   N <- sum(parameters$population)
#
#   # scheduled future doses, by type and dose number
#   events$scheduled_dose <- lapply(X = parameters$N_dose,FUN = function(n_doses){
#     replicate(n = n_doses,expr = individual::TargetedEvent$new(N),simplify = FALSE)
#   })
#   events$scheduled_dose <- setNames(object = events$scheduled_dose,nm = parameters$vaxx_types)
#
#   return(events)
# }

# create_vaccination_dose_listener <- function(variables, parameters, type, dose) {
#   stopifnot(type %in% parameters$vaxx_types)
#   stopifnot( parameters$N_dose[which(parameters$vaxx_types == type)] >= dose )
#   stopifnot( dose > 0 )
#   function (timestep, target) {
#     variables$dose_num$queue_update(values = dose, index = target)
#     variables$dose_time[[dose]]$queue_update(values = timestep, index = target)
#     variables$dose_type[[dose]]$queue_update(values = type, index = target)
#   }
# }
