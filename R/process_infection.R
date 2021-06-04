# --------------------------------------------------
#   infection process for squire transmission model
#   Sean L. Wu (slwood89@gmail.com)
#   June 2021
# --------------------------------------------------


#' @title Infection process (squire model)
#'
#' @description S -> E if infected.
#'
#' @param parameters Model parameters
#' @param variables Model variable
#' @param events Model events
#' @param dt the time step
#' @export
infection_process <- function(parameters, variables, events, dt) {

    return(

      # process without vaccination
      function(timestep) {

        infectious <- variables$states$get_index_of(c("IMild", "IAsymp", "ICase"))

        if (infectious$size() > 0) {

          # Group infection by age
          ages <- variables$discrete_age$get_values(infectious)
          inf_ages <- tabulate(ages, nbins = parameters$N_age)

          # calculate FoI for each age group
          m <- get_contact_matrix(parameters)
          # lambda <- parameters$beta_set[ceiling(timestep * dt)] * rowSums(m %*% diag(inf_ages))
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
