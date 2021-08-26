# --------------------------------------------------------------------------------
#   infection process for vaccination model (multiple doses, no types)
#   Sean L. Wu (slwood89@gmail.com)
#   July 2021
# --------------------------------------------------------------------------------

#' @title Infection process for vaccine model (multi-dose, no types)
#'
#' @description This samples infection events in the susceptible population.
#'
#' @param parameters Model parameters
#' @param variables Model variable
#' @param events Model events
#' @param dt the time step
#' @export
infection_process_vaccine <- function(parameters, variables, events, dt) {

  stopifnot(all(c("states","discrete_age") %in% names(variables)))

  return(

    # process without vaccination
    function(timestep) {

      infectious <- variables$states$get_index_of(c("IMild", "IAsymp", "ICase"))

      if (infectious$size() > 0) {

        # Group infection by age
        ages <- variables$discrete_age$get_values(infectious)
        inf_ages <- tab_bins(a = ages, nbins = parameters$N_age)

        # calculate FoI for each age group
        m <- get_contact_matrix(parameters)
        lambda <- parameters$beta_set[ceiling(timestep * dt)] * as.vector(m %*% inf_ages)

        # Transition from S to E
        susceptible <- variables$states$get_index_of("S")

        # get infection modifier and ages
        ab_titre <- variables$ab_titre$get_values(susceptible)
        infection_efficacy <- vaccine_efficacy_infection(ab_titre = ab_titre,parameters = parameters)
        ages <- variables$discrete_age$get_values(susceptible)

        # FoI for each susceptible based on their age group
        lambda <- lambda[ages]

        # sample infections; individual FoI adjusted by vaccine efficacy
        susceptible$sample(rate = pexp(q = lambda * infection_efficacy * dt))

        # newly infecteds queue the exposure event
        if (susceptible$size() > 0) {
          events$exposure$schedule(susceptible, delay = 0)
        }

      }
    }

  )
}


#' @title C++ infection process for vaccine model (multi-dose, no types)
#'
#' @description This samples infection events in the susceptible population.
#' Calls \code{\link{infection_process_vaccine_cpp_internal}} to return an external pointer object.
#'
#' @param parameters Model parameters
#' @param variables Model variable
#' @param events Model events
#' @param dt the time step
#' @export
infection_process_vaccine_cpp <- function(parameters, variables, events, dt) {

  stopifnot(all(c("states","discrete_age") %in% names(variables)))
  stopifnot("exposure" %in% names(events))

  return(
    infection_process_vaccine_cpp_internal(
      parameters = parameters,
      states = variables$states$.variable,
      discrete_age = variables$discrete_age$.variable,
      ab_titre = variables$ab_titre$.variable,
      exposure = events$exposure$.event,
      dt = dt
    )
  )
}
