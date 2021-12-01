# --------------------------------------------------
#   Translate Ab titre into efficacy
#   and track Ab titre over time
#   Sean L. Wu (slwood89@gmail.com)
#   July 2021
# --------------------------------------------------


#' @title Process that updates the antibody (Ab) titre each time step
#' @description The values in \code{ab_titre} are calculated on the log scale.
#' This process will not calculate decay correctly for `dt > 1` so that is disallowed.
#' @param parameters a list of model parameters
#' @param variables a list of model variables
#' @param vfr an optional vector from [safir::variant_fold_reduction_vector]
#' @param dt time step size
#' @export
vaccine_ab_titre_process <- function(parameters, variables, vfr = NULL, dt) {

  if (!is.null(vfr)) {
    stopifnot(length(vfr) == parameters$time_period / dt)
  }

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

        # if we are supplied with an additional vector for variant fold reduction
        if (!is.null(vfr)) {
          if (vfr[timestep] != 1) {
            # VFR is on linear scale
            new_ab_titre <- exp(new_ab_titre) / vfr[timestep]
            new_ab_titre <- log(new_ab_titre)
          }
        }

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


#' @title Compute vaccine efficacy against infection from Ab titre
#' @param ab_titre a vector of Ab titres
#' @param parameters model parameters.
#' @return a numeric vector, 0 is maximally proective, 1 is maximally unprotective
#' @export
vaccine_efficacy_infection <- function(ab_titre, parameters) {
  # null value is 1
  ef_infection <- rep(1, length(ab_titre))
  if (any(is.finite(ab_titre))) {
    # if some vaccinated individuals with ab titre, calc efficacy for them
    finite_ab <- which(is.finite(ab_titre))
    nt <- exp(ab_titre[finite_ab])
    ef_infection[finite_ab] <- 1 / (1 + exp(-parameters$k * (log10(nt) - log10(parameters$ab_50)))) # reported efficacy in trials
    ef_infection[finite_ab] <- 1 - ef_infection[finite_ab]
    return(ef_infection)
  } else {
    return(ef_infection)
  }
}


#' @title Compute vaccine efficacy against severe disease from Ab titre
#' @description This needs the efficacy against infection because efficacy against severe disease,
#' conditional on breakthrough infection is what safir needs, which is computed as  1 - ((1 - efficacy_disease)/(1 - efficacy_infection)).
#' @param ab_titre a vector of Ab titres
#' @param ef_infection a vector of efficacy against infection from \code{\link{vaccine_efficacy_infection}}
#' @param parameters model parameters
#' @return a numeric vector, 0 is maximally proective, 1 is maximally unprotective
#' @export
vaccine_efficacy_severe <- function(ab_titre, ef_infection, parameters) {
  # null value is 1
  ef_severe <- rep(1, length(ab_titre))
  if (any(is.finite(ab_titre))) {
    # input is on "hazard reduction" scale, convert to efficacy
    # ef_infection <- 1 - ef_infection
    # if some vaccinated individuals with ab titre, calc efficacy for them
    finite_ab <- which(is.finite(ab_titre))
    nt <- exp(ab_titre[finite_ab])
    ef_severe_uncond <- 1 / (1 + exp(-parameters$k * (log10(nt) - log10(parameters$ab_50_severe))))
    ef_severe[finite_ab] <-  1 - ((1 - ef_severe_uncond)/ef_infection[finite_ab]) # 1 - (1 - ef_infection) goes from hazard reduction scale to efficacy, simplifies to ef_infection
    ef_severe[finite_ab] <- 1 - ef_severe[finite_ab]
    return(ef_severe)
  } else {
    return(ef_severe)
  }
}
