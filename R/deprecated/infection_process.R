#' @title Calculating the FOI
#'
#' @description calculating the FOI and infection process
#'
#' @param paramaters Model parameters
#' @param variables Model variable
#' @param events Model events
#' @noRd
infection_process <- function(parameters, variables, events) {

  function(timestep) {
    inf_states <- variables$states$get_index_of(c("IMild", "IAsymp", "ICase"))

    # If IMild = ICase = 0, FOI = 0, i.e. no infected individuals
    if (inf_states$size() > 0) {

      # Group infection by age
      ages <- variables$discrete_age$get_values(inf_states)
      inf_ages <- tabulate(ages, nbins = parameters$N_age)

      # Calculate FoI and use to create probability for each age group
      m <- get_contact_matrix(parameters)

      lambda <- parameters$beta[timestep] * rowSums(m %*% diag(inf_ages))

      # Transition from S to E
      susceptible <- variables$states$get_index_of("S")
      ages <- variables$discrete_age$get_values(susceptible)

      # FOI for each susceptible person
      lambda <- lambda[ages]

      prob_infection  <- 1 - exp(-lambda)

      # infected
      susceptible$sample(rate = prob_infection)

      # newly infecteds queue the exposure event
      if (susceptible$size() > 0) {
        events$exposure$schedule(susceptible, delay = 0)
      }
    }
  }
}
