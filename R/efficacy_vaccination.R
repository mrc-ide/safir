# --------------------------------------------------
#   Translate Ab titre into efficacy
#   and track Ab titre over time
#   Sean L. Wu (slwood89@gmail.com)
#   July 2021
# --------------------------------------------------


#' @title Process that updates the antibody (Ab) titre each time step
#' @description The values in \code{ab_titre} are calculated on the log scale.
#' @param parameters a list of model parameters
#' @param variables a list of model variables
#' @param events a list of model events
#' @param dt time step size
#' @export
vaccine_ab_titre_process <- function(parameters, variables, events, dt) {

  return(
    function(timestep) {

      # only those with at least 1 dose will need to have titre calculated
      vaccinated <- variables$dose_num$get_index_of(set = 0)
      vaccinated <- vaccinated$not()

      if (vaccinated$size() > 0) {

        # for each person we need to know the time since their last dose
        time_since_last_dose <- get_time_since_last_dose(
          timestep = timestep,dt = dt,vaccinated = vaccinated,dose_num = variables$dose_num,dose_time = variables$dose_time,N_phase = parameters$N_phase
        )

        # ceiling to go to next integer day and do not exceed the end of decay rate vector
        time_since_last_dose <- ceiling(time_since_last_dose)
        time_since_last_dose[time_since_last_dose > length(parameters$dr_vec)] <- length(parameters$dr_vec)

        # current Ab titre
        current_ab_titre <- variables$ab_titre$get_values(index = vaccinated)

        # new Ab titre
        new_ab_titre <- current_ab_titre + parameters$dr_vec[time_since_last_dose]

        # schedule an update
        variables$ab_titre$queue_update(values = new_ab_titre, index = vaccinated)

        # vaccine efficacy
        ef_infection <- vaccine_efficacy_infection(ab_titre = current_ab_titre,parameters = parameters)
        ef_severe <- vaccine_efficacy_severe(ab_titre = current_ab_titre,ef_infection = ef_infection,parameters = parameters)

        variables$ef_infection$queue_update(values = ef_infection, index = vaccinated)
        variables$ef_severe$queue_update(values = ef_severe, index = vaccinated)

      }

    }
  )

}


#' @title Calculate the time elapsed in days since each person's last dose
#' @param timestep current time step
#' @param dt size of time step
#' @param vaccinated bitset of vaccinated persons
#' @param dose_num \code{\link[individual]{IntegerVariable}}
#' @param dose_time a list of \code{\link[individual]{IntegerVariable}} objects
#' @param N_phase total number of doses
#' @export
get_time_since_last_dose <- function(timestep, dt, vaccinated, dose_num, dose_time, N_phase) {

  # which dose everybody is on
  vaccinated_dose_num <- dose_num$get_values(vaccinated)

  # output just for the vaccinated persons
  times <- rep(NaN, vaccinated$size())

  for (d in 1:N_phase) {
    times[which(vaccinated_dose_num == d)] <- timestep - dose_time[[d]]$get_values(index = dose_num$get_index_of(set = d))
  }

  times <- times * dt

  return(times)
}


#' @title Compute vaccine efficacy against infection from Ab titre
#' @param ab_titre a vector of Ab titres
#' @param parameters model parameters
#' @export
vaccine_efficacy_infection <- function(ab_titre, parameters) {
  nt <- exp(ab_titre)
  ef_infection <- 1 / (1 + exp(-parameters$k * (log10(nt) - log10(parameters$ab_50)))) # reported efficacy in trials
  ef_infection <- 1 - ef_infection
  return(ef_infection)
}


#' @title Compute vaccine efficacy against severe disease from Ab titre
#' @description This needs the efficacy against infection because efficacy against severe disease,
#' conditional on breakthrough infection is what safir needs, which is computed as  1 - ((1 - efficacy_disease)/(1 - efficacy_infection)).
#' @param ab_titre a vector of Ab titres
#' @param ef_infection a vector of efficacy against infection from \code{\link{vaccine_efficacy_infection}}
#' @param parameters model parameters
#' @export
vaccine_efficacy_severe <- function(ab_titre, ef_infection, parameters) {
  nt <- exp(ab_titre)
  ef_severe_uncond <- 1 / (1 + exp(-parameters$k * (log10(nt) - log10(parameters$ab_50_severe))))
  ef_severe <-  1 - ((1 - ef_severe_uncond)/(1 - ef_infection))
  ef_severe <- 1 - ef_severe
  return(ef_severe)
}
