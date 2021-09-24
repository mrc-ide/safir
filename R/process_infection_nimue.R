# --------------------------------------------------------------------------------
#   infection process for nimue style vaccination model
#   Sean L. Wu (slwood89@gmail.com)
#   June 2021
# --------------------------------------------------------------------------------

#' @title Infection process (nimue vaccine model)
#'
#' @description This samples infection events in the susceptible population. This incorporates
#' the slightly more complex force of infection calculation from the nimue model.
#'
#' @param parameters Model parameters
#' @param variables Model variable
#' @param events Model events
#' @param dt the time step
#' @export
infection_process_nimue <- function(parameters, variables, events, dt) {

  stopifnot(all(c("states","vaccine_states","discrete_age") %in% names(variables)))

  return(

    function(timestep) {

      # current day of simulation
      day <- ceiling(timestep * dt)

      # FoI from contact outside the population
      lambda_external <- parameters$lambda_external[day]

      # infectious classes
      infectious <- variables$states$get_index_of(c("IMild", "IAsymp", "ICase"))

      if (infectious$size() > 0 | lambda_external > 0) {

        susceptible <- variables$states$get_index_of("S")

        # infection by vaccine status
        inf_vaxx <- variables$vaccine_states$get_values(index = infectious)

        # infection by age
        ages <- variables$discrete_age$get_values(infectious)

        # compute cross tab for relative infectiousness, multiply by that matrix, and sum it out
        inf_age_vax <- cross_tab_margins(a = ages,b = inf_vaxx,a_margin = 17,b_margin = 4)
        inf_ages <- rowSums(inf_age_vax * parameters$rel_infectiousness_vaccinated)

        # calculate FoI for each age group
        m <- get_contact_matrix(parameters)
        lambda <- parameters$beta_set[day] * as.vector(m %*% (inf_ages * parameters$rel_infectiousness))

        # FoI for each susceptible person
        sus_vaxx <- variables$vaccine_states$get_values(index = susceptible)
        ages <- variables$discrete_age$get_values(susceptible)

        submat <- matrix(data = NA,nrow = susceptible$size(),ncol = 3)
        submat[, 1] <- day
        submat[, 2] <- ages
        submat[, 3] <- sus_vaxx

        # sample infections
        foi <- (lambda[ages] * parameters$vaccine_efficacy_infection[submat]) + lambda_external
        susceptible$sample(rate = pexp(q = foi * dt))

        # newly infected persons queue the exposure event
        if (susceptible$size() > 0) {
          events$exposure$schedule(susceptible, delay = 0)
        }

      }

    } # end function

  ) # return
}

#' @title C++ infection process (nimue vaccine model)
#'
#' @description Simulates the infection process for the nimue vaccine model.
#' Calls \code{\link{infection_process_nimue_cpp_internal}} to return an external pointer object.
#'
#' @param parameters Model parameters
#' @param variables Model variable
#' @param events Model events
#' @param dt the time step
#' @export
infection_process_nimue_cpp <- function(parameters, variables, events, dt) {

  stopifnot(all(c("states","discrete_age","vaccine_states") %in% names(variables)))
  stopifnot("exposure" %in% names(events))

  return(
    safir::infection_process_nimue_cpp_internal(
      parameters = parameters,
      states = variables$states$.variable,
      vaccine_states = variables$vaccine_states$.variable,
      discrete_age = variables$discrete_age$.variable,
      exposure = events$exposure$.event,
      dt = dt
    )
  )
}

