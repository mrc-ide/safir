#' @title Process that updates the antibody (Ab) titre each time step
#' @description The values in `ab_titre` are calculated on the log scale.
#' This process will not calculate decay correctly for `dt > 1` so that is disallowed.
#' @param parameters a list of model parameters
#' @param variables a list of model variables
#' @param dt time step size
#' @export
vaccine_ab_titre_process <- function(parameters, variables, dt) {

  stopifnot(dt <= 1)

  return(
    function(timestep) {

      # only those with at least 1 dose will need to have titre calculated
      vaccinated <- variables$dose_num$get_index_of(set = 0)
      vaccinated$not(inplace = TRUE)

      if (vaccinated$size() > 0) {

        # for each person we need to know the time since their last dose
        time_since_last_dose <- get_time_since_last_dose(timestep = timestep, dt = dt, vaccinated = vaccinated, dose_time = variables$dose_time)

        # ceiling to go to next integer day and do not exceed the end of decay rate vector
        time_since_last_dose <- ceiling(time_since_last_dose)
        time_since_last_dose[time_since_last_dose > length(parameters$dr_vec)] <- length(parameters$dr_vec)

        # current Ab titre
        current_ab_titre <- variables$ab_titre$get_values(index = vaccinated)

        # new Ab titre
        new_ab_titre <- current_ab_titre + (parameters$dr_vec[time_since_last_dose] * dt)

        # schedule an update
        variables$ab_titre$queue_update(values = new_ab_titre, index = vaccinated)

      }

    }
  )

}

#' @title Calculate the time elapsed in days since each person's last dose
#' @param timestep current time step
#' @param dt size of time step
#' @param vaccinated bitset of vaccinated persons
#' @param dose_time a list of \code{\link[individual]{IntegerVariable}} objects
#' @export
get_time_since_last_dose <- function(timestep, dt, vaccinated, dose_time) {

  # output just for the vaccinated persons
  times <- rep(NaN, vaccinated$size())
  times <- timestep - dose_time$get_values(index = vaccinated)
  times <- times * dt

  return(times)
}



#' @title Process that updates the antibody (Ab) titre each time step for model vaccine model with natural immunity
#' @description The values in `ab_titre` are calculated on the log scale.
#' This process will not calculate decay correctly for `dt > 1` so that is disallowed.
#' @param parameters a list of model parameters
#' @param variables a list of model variables
#' @param dt time step size
#' @export
natural_immunity_ab_titre_process <- function(parameters, variables, dt) {

  return(
    function(timestep) {

      # persons who have been vaccinated or infected will have antibody titre
      vaccinated <- variables$dose_num$get_index_of(set = 0)$not(inplace = TRUE)
      infected <- variables$inf_num$get_index_of(set = 0)$not(inplace = TRUE)

      vaccinated_or_infected <- vaccinated$or(infected)

      if (vaccinated_or_infected$size() > 0) {

        # for each person we need to know the time since their last dose
        time_since_last_dose_or_infection <- get_time_since_last_dose_or_infection(
          timestep = timestep, dt = dt, vaccinated_or_infected = vaccinated_or_infected,
          dose_time = variables$dose_time, inf_time = variables$inf_time
        )

        # ceiling to go to next integer day and do not exceed the end of decay rate vector
        time_since_last_dose_or_infection <- ceiling(time_since_last_dose_or_infection)
        time_since_last_dose_or_infection[time_since_last_dose_or_infection > length(parameters$dr_vec)] <- length(parameters$dr_vec)

        # current Ab titre
        current_ab_titre <- variables$ab_titre$get_values(index = vaccinated_or_infected)

        # new Ab titre
        new_ab_titre <- current_ab_titre + (parameters$dr_vec[time_since_last_dose_or_infection] * dt)

        # schedule an update
        variables$ab_titre$queue_update(values = new_ab_titre, index = vaccinated_or_infected)

      }

    }
  )

}

#' @title Calculate the time elapsed in days since each person's last dose
#' @param timestep current time step
#' @param dt size of time step
#' @param vaccinated_or_infected [individual::Bitset] of vaccinated persons
#' @param dose_time an [individual::IntegerVariable] object
#' @param inf_time an [individual::IntegerVariable] object
#' @export
get_time_since_last_dose_or_infection <- function(timestep, dt, vaccinated_or_infected, dose_time, inf_time) {

  last_dose_tt <- dose_time$get_values(vaccinated_or_infected)
  last_inf_tt <- inf_time$get_values(vaccinated_or_infected)

  # the max of last dose and last inf time will give the time of last ab boost
  times <- timestep - pmax(last_dose_tt, last_inf_tt)
  times <- times * dt

  return(times)
}
