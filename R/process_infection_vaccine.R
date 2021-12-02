# --------------------------------------------------------------------------------
#   infection process for vaccination model (multiple doses, no types)
#   Sean L. Wu (slwood89@gmail.com)
#   July 2021
# --------------------------------------------------------------------------------

#' @title Infection process for vaccine model (multi-dose, no types)
#'
#' @description This samples infection events in the susceptible population.
#'
#' @param parameters Model parameters
#' @param variables Model variable
#' @param events Model events
#' @param dt the time step
#' @export
infection_process_vaccine <- function(parameters, variables, events, dt) {

  stopifnot(all(c("states","discrete_age") %in% names(variables)))

  if (parameters$nt_efficacy_transmission) {
    get_inf_ages <- function(infection_bset, variables, parameters, timestep) {
      ages <- variables$discrete_age$get_values(infection_bset)
      nat_values <- variables$ab_titre$get_values(infection_bset)
      inf_wt <- vaccine_efficacy_transmission_cpp(ab_titre = nat_values, parameters = parameters, timestep = timestep)
      inf_ages <- tab_bins_weighted(a = ages, wt = inf_wt,  nbins = parameters$N_age)
      return(inf_ages)
    }
  } else {
    get_inf_ages <- function(infection_bset, variables, parameters, timestep) {
      ages <- variables$discrete_age$get_values(infection_bset)
      inf_ages <- tab_bins(a = ages, nbins = parameters$N_age)
      return(inf_ages)
    }
  }

  return(

    # process without vaccination
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

          # group infectious persons by age
          inf_ages <- get_inf_ages(infection_bset = infectious, variables = variables, parameters = parameters, timestep = day)

          # calculate FoI on each susceptible age group
          m <- get_contact_matrix(parameters)
          lambda_age <- parameters$beta_set[day] * as.vector(m %*% inf_ages)

          # get infection modifier and ages
          ab_titre <- variables$ab_titre$get_values(susceptible)
          infection_efficacy <- vaccine_efficacy_infection_cpp(ab_titre = ab_titre,parameters = parameters, timestep = timestep)
          ages <- variables$discrete_age$get_values(susceptible)

          # FoI for each susceptible based on their age group
          lambda <- lambda + (lambda_age[ages] * infection_efficacy)

        }

        # sample infection events in susceptible population
        susceptible$sample(rate = pexp(q = lambda * dt))

        # queue the exposure event
        events$exposure$schedule(susceptible, delay = 0)

      } # end if S > 0

    } # end process fn

  ) # end return

}


#' @title C++ infection process for vaccine model (multi-dose, no types)
#'
#' @description This samples infection events in the susceptible population.
#' Calls \code{\link{infection_process_vaccine_cpp_internal}} to return an external pointer object.
#'
#' @param parameters Model parameters
#' @param variables Model variable
#' @param events Model events
#' @param dt the time step
#' @export
infection_process_vaccine_cpp <- function(parameters, variables, events, dt) {

  stopifnot(all(c("states","discrete_age") %in% names(variables)))
  stopifnot("exposure" %in% names(events))

  return(
    infection_process_vaccine_cpp_internal(
      parameters = parameters,
      states = variables$states$.variable,
      discrete_age = variables$discrete_age$.variable,
      ab_titre = variables$ab_titre$.variable,
      exposure = events$exposure$.event,
      dt = dt
    )
  )
}
