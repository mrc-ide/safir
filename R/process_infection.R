# --------------------------------------------------
#   infection process for squire transmission model
#   Sean L. Wu (slwood89@gmail.com)
#   June 2021
# --------------------------------------------------


#' @title Infection process (squire transmission model)
#'
#' @description simulates the process of infection.
#'
#' @param parameters Model parameters
#' @param variables a list of model variables, the output of [create_variables]
#' @param events a list of [individual::TargetedEvent], the output of [create_events]
#' @param dt the time step
#' @export
infection_process <- function(parameters, variables, events, dt) {

    return(

      # process without vaccination
      function(timestep) {

        # current day of simulation
        day <- ceiling(timestep * dt)

        # FoI from contact outside the population
        lambda_external <- parameters$lambda_external[day]

        # infectious classes
        infectious <- variables$states$get_index_of(c("IMild", "IAsymp", "ICase"))

        if (infectious$size() > 0 | lambda_external > 0) {

          # Group infection by age
          ages <- variables$discrete_age$get_values(infectious)
          inf_ages <- tab_bins(a = ages, nbins = parameters$N_age)

          # calculate FoI for each age group
          m <- get_contact_matrix(parameters)
          lambda <- parameters$beta_set[day] * as.vector(m %*% inf_ages)

          # Transition from S to E
          susceptible <- variables$states$get_index_of("S")
          ages <- variables$discrete_age$get_values(susceptible)

          # FoI for each susceptible person
          lambda <- lambda[ages]

          # infected
          susceptible$sample(rate = pexp(q = (lambda + lambda_external) * dt))

          # newly infecteds queue the exposure event
          if (susceptible$size() > 0) {
            events$exposure$schedule(susceptible, delay = 0)
          }

        }
      }

    )
}


#' @title C++ infection process (squire transmission model)
#'
#' @description Simulates the infection process for the squire transmission model.
#' Calls \code{\link{infection_process_cpp_internal}} to return an external pointer object.
#'
#' @param parameters Model parameters
#' @param variables a list of model variables, the output of [create_variables]
#' @param events a list of [individual::TargetedEvent], the output of [create_events]
#' @param dt the time step
#' @export
infection_process_cpp <- function(parameters, variables, events, dt) {

  stopifnot(all(c("states","discrete_age") %in% names(variables)))
  stopifnot("exposure" %in% names(events))

  return(
    infection_process_cpp_internal(
      parameters = parameters,
      states = variables$states$.variable,
      discrete_age = variables$discrete_age$.variable,
      exposure = events$exposure$.event,
      dt = dt
    )
  )
}
