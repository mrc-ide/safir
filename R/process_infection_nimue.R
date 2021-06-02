#' @title Transitions from S (nimue vaccine model)
#'
#' @description S -> E if infected; can also become vaccinated. This incorporates
#' the slightly more complex force of infection calculation from the nimue model.
#' This function assumes an individual can become infected and vaccinated in a
#' time step.
#'
#' @param parameters Model parameters
#' @param variables Model variable
#' @param events Model events
#' @param dt the time step
#' @export
infection_process_nimue <- function(parameters, variables, events, dt) {

  return(

    function(timestep) {

      # infection
      infectious <- variables$states$get_index_of(c("IMild", "IAsymp", "ICase"))
      if (infectious$size() > 0) {

        # Group infection by age
        ages <- variables$discrete_age$get_values(infectious)
        inf_ages <- tabulate(ages, nbins = parameters$N_age)

        # calculate FoI for each age group
        m <- get_contact_matrix(parameters)
        lambda <- parameters$beta_set[ceiling(timestep * dt)] * rowSums(m %*% diag(inf_ages))

        # Transition from S to E
        susceptible <- variables$states$get_index_of("S")
        ages <- variables$discrete_age$get_values(susceptible)

        # FoI for each susceptible person
        lambda <- lambda[ages]
        susceptible$sample(rate = pexp(q = lambda * dt))

        # newly infecteds queue the exposure event
        if (susceptible$size() > 0) {
          events$exposure$schedule(susceptible, delay = 0)
        }

      }

      # vaccination
      if (variables$eligible$size() > 0) {

        vaccinated <- variables$eligible$copy()
        vaccinated$sample(rate = pexp(q = variables$vr * dt))

        if (vaccinated$size() > 0) {
          events$v0_to_v1v2$schedule(vaccinated, delay = 0)
        }

      }


    } # end function

  ) # return
}
