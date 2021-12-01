#' @title Compute vaccine efficacy against infection from Ab titre
#' @param ab_titre a vector of Ab titres
#' @param parameters model parameters.
#' @param timestep the time step
#' @return a numeric vector, 0 is maximally proective, 1 is maximally unprotective
#' @export
vaccine_efficacy_infection <- function(ab_titre, parameters, timestep) {

  if (is.null(parameters$vfr)) {
    vfr <- 1
  } else {
    vfr <- parameters$vfr[timestep]
  }

  # null value is 1
  ef_infection <- rep(1, length(ab_titre))

  if (any(is.finite(ab_titre))) {
    # if some vaccinated individuals with ab titre, calc efficacy for them
    finite_ab <- which(is.finite(ab_titre))
    nt <- exp(ab_titre[finite_ab])
    nt <- pmax(.Machine$double.eps, nt / vfr)
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
#' @param timestep the time step
#' @return a numeric vector, 0 is maximally proective, 1 is maximally unprotective
#' @export
vaccine_efficacy_severe <- function(ab_titre, ef_infection, parameters, timestep) {

  if (is.null(parameters$vfr)) {
    vfr <- 1
  } else {
    vfr <- parameters$vfr[timestep]
  }

  # null value is 1
  ef_severe <- rep(1, length(ab_titre))

  if (any(is.finite(ab_titre))) {
    # input is on "hazard reduction" scale, convert to efficacy
    # ef_infection <- 1 - ef_infection
    # if some vaccinated individuals with ab titre, calc efficacy for them
    finite_ab <- which(is.finite(ab_titre))
    nt <- exp(ab_titre[finite_ab])
    nt <- pmax(.Machine$double.eps, nt / vfr)
    ef_severe_uncond <- 1 / (1 + exp(-parameters$k * (log10(nt) - log10(parameters$ab_50_severe))))
    ef_severe[finite_ab] <-  1 - ((1 - ef_severe_uncond)/ef_infection[finite_ab]) # 1 - (1 - ef_infection) goes from hazard reduction scale to efficacy, simplifies to ef_infection
    ef_severe[finite_ab] <- 1 - ef_severe[finite_ab]
    return(ef_severe)
  } else {
    return(ef_severe)
  }
}


#' @title Compute vaccine efficacy against onward transmission from Ab titre
#' @param ab_titre a vector of Ab titres
#' @param parameters model parameters.
#' @param timestep the time step
#' @return a numeric vector, 0 is maximally protective, 1 is maximally unprotective
#' @export
vaccine_efficacy_transmission <- function(ab_titre, parameters, timestep) {

  if (is.null(parameters$vfr)) {
    vfr <- 1
  } else {
    vfr <- parameters$vfr[timestep]
  }

  # null value is 1
  ef_transmission <- rep(1, length(ab_titre))

  if (any(is.finite(ab_titre))) {
    # if some vaccinated individuals with ab titre, calc efficacy for them
    finite_ab <- which(is.finite(ab_titre))
    nt <- exp(ab_titre[finite_ab])
    nt <- pmax(.Machine$double.eps, nt / vfr)
    ef_transmission[finite_ab] <- 1 / (1 + exp(-parameters$k * (log10(nt / parameters$nt_transmission_factor) - log10(parameters$ab_50)))) # reported efficacy in trials
    ef_transmission[finite_ab] <- 1 - ef_transmission[finite_ab]
    return(ef_transmission)
  } else {
    return(ef_transmission)
  }
}
