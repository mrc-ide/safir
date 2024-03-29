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

      # susceptible persons
      susceptible <- variables$states$get_index_of("S")

      if (susceptible$size() > 0) {

        # FoI for each susceptible from external contacts
        lambda <- rep(x = lambda_external, times = susceptible$size())

        # FoI contribution from transmission
        if (infectious$size() > 0) {

          # group infectious persons by vaccine status and age
          vaxx <- variables$vaccine_states$get_values(index = infectious)
          ages <- variables$discrete_age$get_values(infectious)

          # compute cross tab for relative infectiousness, multiply by that matrix, and sum it out
          inf_age_vax <- cross_tab_margins(a = ages,b = vaxx,a_margin = 17,b_margin = 4)
          inf_ages <- rowSums(inf_age_vax * parameters$rel_infectiousness_vaccinated)

          # calculate FoI on each susceptible age group
          m <- get_contact_matrix(parameters)
          lambda_age <- parameters$beta_set[day] * as.vector(m %*% (inf_ages * parameters$rel_infectiousness))

          # FoI for each susceptible person from transmission accounting for age and vaccine status
          sus_vaxx <- variables$vaccine_states$get_values(index = susceptible)
          ages <- variables$discrete_age$get_values(susceptible)

          submat <- matrix(data = NA,nrow = susceptible$size(),ncol = 3)
          submat[, 1] <- day
          submat[, 2] <- ages
          submat[, 3] <- sus_vaxx

          # FoI on infectives from transmission
          lambda <- lambda + (lambda_age[ages] * parameters$vaccine_efficacy_infection[submat])

        }

        # sample infection events in susceptible population
        susceptible$sample(rate = pexp(q = lambda * dt))

        # queue the exposure event
        events$exposure$schedule(susceptible, delay = 0)

      } # end if S > 0

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

