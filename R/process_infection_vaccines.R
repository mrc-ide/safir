# --------------------------------------------------------------------------------
#   infection process for vaccination model (multiple doses, no types)
#   Sean L. Wu (slwood89@gmail.com)
#   July 2021
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
infection_process_vaccine <- function(parameters, variables, events, dt) {

  stopifnot(all(c("states","vaccine_states","discrete_age") %in% names(variables)))

  # doses_vec <- as.character(0:parameters$N_phase)
  # doses_bset <- replicate(n = length(doses_vec),expr = {NULL})

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
        ages <- variables$discrete_age$get_values(susceptible)

        # FoI for each susceptible person
        lambda <- lambda[ages]

        # infected
        susceptible$sample(rate = pexp(q = lambda * dt))

        # newly infecteds queue the exposure event
        if (susceptible$size() > 0) {
          events$exposure$schedule(susceptible, delay = 0)
        }

      }
    }

  )
}
