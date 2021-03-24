#' @title Calculating the FOI
#'
#' @description calculating the FOI and infection process
#'
#' @param human the human handle
#' @param states a list of states in the model
#' @param discrete_age variable
#' @param exposure event for covid exposure
#' @param contact_matrix_set contact matrix set
#' @noRd
infection_process <- function(human, states, discrete_age, exposure,
                              contact_matrix_set) {

  function(api) {
    pars <- api$get_parameters()

    inf_states <- api$get_state(human, states$IMild, states$IAsymp,
                                states$ICase)

    # If IMild = ICase = 0, FOI = 0, i.e. no infected individuals
    if (length(inf_states) > 0) {

      # Group infection by age
      ages <- api$get_variable(human, discrete_age, inf_states)
      inf_ages <- tabulate(ages, nbins = pars$N_age)

      # Calculate FoI and use to create probability for each age group
      m <- get_contact_matrix(contact_matrix_set)

      lambda <- pars$beta[api$get_timestep()] * rowSums(m %*% diag(inf_ages))

      # Transition from S to E
      susceptible <- api$get_state(human, states$S)
      ages <- api$get_variable(human, discrete_age, susceptible)

      # FOI for each susceptible person
      lambda <- lambda[ages]
      prob_infection  <- 1 - exp(-lambda)

      # infected
      infected <- bernoulli_multi_p(prob_infection)

      # if infections then
      if(sum(infected) > 0) {
        api$schedule(
          event = exposure,
          target = susceptible[infected],
          delay = 0 # i.e. happens now
        )
      }
    }
  }
}
