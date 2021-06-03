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

      infectious <- variables$states$get_index_of(c("IMild", "IAsymp", "ICase"))

      if (infectious$size() > 0) {

        susceptible <- variables$states$get_index_of("S")

        # infection by vaccine status
        inf_vaxx <- variables$vaccine_states$get_values(index = infectious)

        # infection by age
        ages <- variables$discrete_age$get_values(infectious)

        # compute cross tab for relative infectiousness, multiply by that matrix, and sum it out
        inf_age_vax <- cross_tab_margins(a = ages,b = inf_vaxx,a_margin = 1:17,b_margin = 1:4)
        inf_ages <- rowSums(inf_age_vax * parameters$rel_infectiousness_vaccinated)

        # calculate FoI for each age group
        m <- get_contact_matrix(parameters)
        lambda <- parameters$beta_set[ceiling(timestep * dt)] * rowSums(m %*% diag(inf_ages) %*% diag(parameters$rel_infectiousness))

        # FoI for each susceptible person
        ages <- variables$discrete_age$get_values(susceptible)
        lambda <- lambda[ages]
        susceptible$sample(rate = pexp(q = lambda * dt))

        # newly infecteds queue the exposure event
        if (susceptible$size() > 0) {
          events$exposure$schedule(susceptible, delay = 0)
        }

      }

    } # end function

  ) # return
}
