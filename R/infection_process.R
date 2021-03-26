#' @title Calculating the FOI
#'
#' @description calculating the FOI and infection process
#'
#' @param exposure event for covid exposure
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
      m <- get_contact_matrix(parameters$mix_mat_set)

      lambda <- parameters$beta[timestep] * rowSums(m %*% diag(inf_ages))

      # Transition from S to E
      susceptible <- variables$states$get_index_of("S")
      ages <- variables$discrete_age$get_values(susceptible)

      # FOI for each susceptible person
      lambda <- lambda[ages]

      prob_infection  <- 1 - exp(-lambda)

      # infected
      infected <- bernoulli_multi_p(prob_infection)

      # if infections then
      if(sum(infected) > 0) {
        to_infect <- individual::filter_bitset(susceptible, which(infected))
        events$exposure$schedule(to_infect, delay = 0)
      }
    }
  }
}
