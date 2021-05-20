# --------------------------------------------------
#   infection process
#   May 2021
#   1. infection_process
# --------------------------------------------------


#' @title Transitions from S
#'
#' @description S -> E if infected; can also become vaccinated.
#'
#' @param parameters Model parameters
#' @param variables Model variable
#' @param events Model events
#' @param dt the time step
#' @param vaccines vaccine parameters
#' @export
infection_process <- function(parameters, variables, events, dt, vaccines = NULL) {

  if (is.null(vaccines)) {

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
          lambda <- parameters$beta[ceiling(timestep * dt)] * rowSums(m %*% diag(inf_ages))

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

  } else {
    stop("not implemented yet!")
  }

}
