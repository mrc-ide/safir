# --------------------------------------------------
#   Translate Ab titre into efficacy
#   and track Ab titre over time
#   Sean L. Wu (slwood89@gmail.com)
#   July 2021
# --------------------------------------------------


#' @title Process that updates the antibody (Ab) titre each time step
#' @param ab_titre a \code{\link[individual]{DoubleVariable}} giving the antibody titre for the population
#' @param dose_time a list of \code{\link[individual]{IntegerVariable}} objects giving the time of doses for the population
#' @param dose_num a \code{\link[individual]{IntegerVariable}} giving the dose number of the population
#' @export
vaccine_ab_titre_process <- function(parameters, variables, events, dt) {

  dose_num <- variables$dose_num
  dose_time <- variables$dose_time
  ab_titre <- variables$ab_titre

  return(
    function(timestep) {

      # only those with at least 1 dose will need to have titre calculated
      vaccinated <- dose_num$get_index_of(set = 0)
      vaccinated <- vaccinated$not()

      if (vaccinated$size() > 0) {

        # for each person we need to know the time since their last dose
        time_since_last_dose <- get_time_since_last_dose(
          timestep = timestep,dt = dt,vaccinated = vaccinated,dose_num = variables$dose_num,dose_time = variables$dose_time,N_phase = parameters$N_phase
        )

        # current Ab titre
        current_ab_values <- ab_titre$get_values(index = vaccinated)

        # apply the discrete difference operator
        # new_ab_values <- blah
        new_ab_values <- pmax(x - (x * 0.01 * dt), 0) # fake dynamics

        # schedule an update
        ab_titre$queue_update(values = new_ab_values, index = vaccinated)
      }

    }
  )

}


#' @title Calculate the time elapsed in days since each person's last dose
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
#' @param ab_titre a \code{\link[individual]{DoubleVariable}} containing values of Ab titre
#' @param who a \code{\link[individual]{Bitset}} telling who we should calculate efficacy against infection for
#' @export
vaccine_efficacy_infection <- function(ab_titre, who) {
  ab_values <- ab_titre$get_values(index = who)
  return(exp(-ab_values/3))
}


#' @title Compute vaccine efficacy against severe disease from Ab titre
#' @param ab_titre a \code{\link[individual]{DoubleVariable}} containing values of Ab titre
#' @param who a \code{\link[individual]{Bitset}} telling who we should calculate efficacy against severe infection for
#' @export
vaccine_efficacy_severe <- function(ab_titre, who) {
  ab_values <- ab_titre$get_values(index = who)
  return(exp(-ab_values/2))
}
