#' @title Compute vaccine efficacy against infection from Ab titre
#' @param nat a vector of NAT values
#' @param parameters model parameters.
#' @param day the current day
#' @return a numeric vector, 0 is maximally proective, 1 is maximally unprotective
#' @export
vaccine_efficacy_infection <- function(nat, parameters, day) {

  # null value is 1
  ef_infection <- rep(1, length(nat))

  nonzero_nat <- which(nat > 0)

  ef_infection[nonzero_nat] <- 1 / (1 + exp(-parameters$k * (log10(nat[nonzero_nat]) - log10(parameters$ab_50)))) # reported efficacy in trials
  ef_infection[nonzero_nat] <- 1 - ef_infection[nonzero_nat]
  return(ef_infection)
}


#' @title Compute vaccine efficacy against severe disease from Ab titre
#' @description This needs the efficacy against infection because efficacy against severe disease,
#' conditional on breakthrough infection is what safir needs, which is computed as  1 - ((1 - efficacy_disease)/(1 - efficacy_infection)).
#' @param nat a vector of NAT values
#' @param ef_infection a vector of efficacy against infection from \code{\link{vaccine_efficacy_infection}}
#' @param parameters model parameters
#' @param day the current day
#' @return a numeric vector, 0 is maximally proective, 1 is maximally unprotective
#' @export
vaccine_efficacy_severe <- function(nat, ef_infection, parameters, day) {

  # null value is 1
  ef_severe <- rep(1,length(nat))

  nonzero_nat <- which(nat > 0)

  ef_severe_uncond <- 1 / (1 + exp(-parameters$k * (log10(nat[nonzero_nat]) - log10(parameters$ab_50_severe))))
  ef_severe[nonzero_nat] <-  1 - ((1 - ef_severe_uncond)/ef_infection[nonzero_nat]) # 1 - (1 - ef_infection) goes from hazard reduction scale to efficacy, simplifies to ef_infection
  ef_severe[nonzero_nat] <- 1 - ef_severe[nonzero_nat]
  return(ef_severe)
}


#' @title Compute vaccine efficacy against onward transmission from Ab titre
#' @param nat a vector of NAT values
#' @param parameters model parameters.
#' @param day the current day
#' @return a numeric vector, 0 is maximally protective, 1 is maximally unprotective
#' @export
vaccine_efficacy_transmission <- function(nat, parameters, day) {

  # null value is 1
  ef_transmission <- rep(1, length(nat))

  nonzero_nat <- which(nat > 0)

  ef_transmission[nonzero_nat] <- 1 / (1 + exp(-parameters$k * (log10(nat[nonzero_nat] / parameters$nt_transmission_factor) - log10(parameters$ab_50)))) # reported efficacy in trials
  ef_transmission[nonzero_nat] <- 1 - ef_transmission[nonzero_nat]

  return(ef_transmission)
}


#' @title Make function to calculate population NAT
#' @description The returned function takes arguments `index` (a bitset) and `day` and returns the NAT of persons in `index` on the linear scale.
#' @param variables a named list
#' @param parameters a named list
#' @export
make_calculate_nat <- function(variables, parameters) {
  stopifnot(!is.null(variables$ab_titre))
  stopifnot(inherits(variables$ab_titre, "DoubleVariable"))
  if (is.null(variables$ab_titre_inf)) {

    stop("currently safir only correctly simulates models with seperate tracking of vaccine and infection derived NAT values")

  } else {
    stopifnot(inherits(variables$ab_titre_inf, "DoubleVariable"))

    # for new variant proof VFR modeling
    if (!is.null(parameters$vp_time)) {

      dt <- parameters$dt

      calculate_nat <- function(index, day) {

        # Get VFR parameter, if null assign 1
        if (is.null(parameters$vfr)) {
          vfr <- 1
        } else {
          vfr <- parameters$vfr[day]
        }

        nat_vaccine <- exp(variables$ab_titre$get_values(index))
        nat_infection <- exp(variables$ab_titre_inf$get_values(index))

        # nat_infection should always be scaled by vfr regardless
        nat_infection <-  nat_infection / vfr

        # find individuals who were vaccinated when variant proof vaccine was on
        dose_times <- ceiling(variables$dose_time$get_values(index) * dt)

        # Find those in index that have been vaccinated (dose_num not equal to 0)
        been_vaccinated_vec <- which(variables$dose_num$get_values(index) != 0)

        # who are those that were vaccinated during this time
        vp_index <- which(parameters$vp_time[dose_times[been_vaccinated_vec]] == 1)

        # people who have been vaccinated and whose dose_times are not during variant proof window
        been_vaccinated_with_vp <- been_vaccinated_vec[vp_index]

        # apply vfr to everyone not vaccinated with VP
        vfr_indexes <- setdiff(seq_along(nat_vaccine), been_vaccinated_with_vp)
        nat_vaccine[vfr_indexes] <-  nat_vaccine[vfr_indexes] / vfr

        # and combine these for overall NAT
        nt <- nat_vaccine + nat_infection

        # and now do the pmax comparison
        nt <- pmax(.Machine$double.eps, nt)

        return(nt)
      }

    } else {
      # existing setup

      calculate_nat <- function(index, day) {

        # Get VFR parameter, if null assign 1
        if (is.null(parameters$vfr)) {
          vfr <- 1
        } else {
          vfr <- parameters$vfr[day]
        }

        nat_vaccine <- exp(variables$ab_titre$get_values(index))
        nat_infection <- exp(variables$ab_titre_inf$get_values(index))

        nt <- nat_vaccine + nat_infection
        nt <- pmax(.Machine$double.eps, nt / vfr)
        return(nt)
      }

    }

  }
  return(calculate_nat)
}
