#' @title Compute vaccine efficacy against infection from Ab titre
#' @param variables a list of variables that include the ab_titres
#' @param parameters model parameters.
#' @param day the current day
#' @return a numeric vector, 0 is maximally proective, 1 is maximally unprotective
#' @export
vaccine_efficacy_infection <- function(variables, parameters, day) {

  # Check if the ab_titre vectors are not null
  stopifnot(!is.null(variables$ab_titre))
  stopifnot(!is.null(variables$ab_titre_inf))
  stopifnot(inherits(variables$ab_titre, "DoubleVariable"))

  # Get VFR parameter, if null assign 1
  if (is.null(parameters$vfr)) {
    vfr <- 1
  } else {
    vfr <- parameters$vfr[day]
  }

  # null value is 1
  ef_infection <- rep(1,variables$ab_titre)


  if (any(is.finite(variables$ab_titre))){

    # if some vaccinated individuals with ab titre, calc efficacy for them
    finite_ab <- which(variables$ab_titre)

    nat_vaccine <- exp(variables$ab_titre$get_values(index))
    nat_infection <- exp(variables$ab_titre_inf$get_values(index))

    vp_on <- parameters$vp_time[day]

    if(vp_on == 0 ){  # If there is not a variant proof vaccine

      nt <- nat_vaccine + nat_infection
      nt <- pmax(.Machine$double.eps, nt / vfr)

    } else if (vp_on ==1 ){

      nat_infection <-  pmax(.Machine$double.eps, nat_infection / vfr)  # Only infection induced antibodies are scaled down by vfr
      nt <- nat_vaccine + nat_infection
    }

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
#' @param variables a list of variables that include the ab_titres
#' @param ef_infection a vector of efficacy against infection from \code{\link{vaccine_efficacy_infection}}
#' @param parameters model parameters
#' @param day the current day
#' @return a numeric vector, 0 is maximally proective, 1 is maximally unprotective
#' @export
vaccine_efficacy_severe <- function(ab_titre, ef_infection, parameters, day) {


  # Check if the ab_titre vectors are not null
  stopifnot(!is.null(variables$ab_titre))
  stopifnot(!is.null(variables$ab_titre_inf))
  stopifnot(inherits(variables$ab_titre, "DoubleVariable"))

  # Get VFR parameter, if null assign 1
  if (is.null(parameters$vfr)) {
    vfr <- 1
  } else {
    vfr <- parameters$vfr[day]
  }

  # null value is 1
  ef_severe <- rep(1, length(variables$ab_titre))

  if (any(is.finite(variables$ab_titre))){

    # if some vaccinated individuals with ab titre, calc efficacy for them
    finite_ab <- which(variables$ab_titre)

    nat_vaccine <- exp(variables$ab_titre$get_values(index))
    nat_infection <- exp(variables$ab_titre_inf$get_values(index))

    vp_on <- parameters$vp_time[day]

    if(vp_on == 0 ){  # If there is not a variant proof vaccine

      nt <- nat_vaccine + nat_infection
      nt <- pmax(.Machine$double.eps, nt / vfr)

    } else if (vp_on ==1 ){

      nat_infection <-  pmax(.Machine$double.eps, nat_infection / vfr)  # Only infection induced antibodies are scaled down by vfr
      nt <- nat_vaccine + nat_infection
    }


    ef_severe_uncond <- 1 / (1 + exp(-parameters$k * (log10(nt) - log10(parameters$ab_50_severe))))
    ef_severe[finite_ab] <-  1 - ((1 - ef_severe_uncond)/ef_infection[finite_ab]) # 1 - (1 - ef_infection) goes from hazard reduction scale to efficacy, simplifies to ef_infection
    ef_severe[finite_ab] <- 1 - ef_severe[finite_ab]
    return(ef_severe)
  } else {
    return(ef_severe)
  }
}


#' @title Compute vaccine efficacy against onward transmission from Ab titre
#' @param variables a list of variables that include the ab_titres
#' @param parameters model parameters.
#' @param day the current day
#' @return a numeric vector, 0 is maximally protective, 1 is maximally unprotective
#' @export
vaccine_efficacy_transmission <- function(ab_titre, parameters, day) {

  # Check if the ab_titre vectors are not null
  stopifnot(!is.null(variables$ab_titre))
  stopifnot(!is.null(variables$ab_titre_inf))
  stopifnot(inherits(variables$ab_titre, "DoubleVariable"))

  # Get VFR parameter, if null assign 1
  if (is.null(parameters$vfr)) {
    vfr <- 1
  } else {
    vfr <- parameters$vfr[day]
  }

  # null value is 1
  ef_transmission <- rep(1, length(variables$ab_titre))

  if (any(is.finite(variables$ab_titre))){
    # if some vaccinated individuals with ab titre, calc efficacy for them
    finite_ab <- which(variables$ab_titre)

    nat_vaccine <- exp(variables$ab_titre$get_values(index))
    nat_infection <- exp(variables$ab_titre_inf$get_values(index))

    vp_on <- parameters$vp_time[day]

    if(vp_on == 0 ){  # If there is not a variant proof vaccine

      nt <- nat_vaccine + nat_infection
      nt <- pmax(.Machine$double.eps, nt / vfr)

    } else if (vp_on ==1 ){

      nat_infection <-  pmax(.Machine$double.eps, nat_infection / vfr)  # Only infection induced antibodies are scaled down by vfr
      nt <- nat_vaccine + nat_infection
    }


    ef_transmission[finite_ab] <- 1 / (1 + exp(-parameters$k * (log10(nt / parameters$nt_transmission_factor) - log10(parameters$ab_50)))) # reported efficacy in trials
    ef_transmission[finite_ab] <- 1 - ef_transmission[finite_ab]
    return(ef_transmission)
  } else {
    return(ef_transmission)
  }
}

