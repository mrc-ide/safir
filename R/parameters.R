# --------------------------------------------------
#   model parameters
#   May 2021
#   1. get_parameters
#   2. get_population
#   3. get_country
#   4. get_asymptomatic
#   5. interp_vector_constant
#   6. interp_matrix_list_constant
#   7. interp_input_par
# --------------------------------------------------


#' @title Get parameters from squire model
#'
#' @inheritParams squire::parameters_explicit_SEEIR
#' @param iso3c Character for country iso3c
#' @param max_age the maximum age for humans
#' @param ... Other parameters for [squire::parameters_explicit_SEEIR]
#'
#' @return squire model parameters
#' @export
get_parameters <- function(iso3c = NULL,
                           population = NULL,
                           contact_matrix_set = NULL,
                           time_period = 365,
                           max_age = 100,
                           ...) {

  # if missing a contact matrix but have iso3c then use that
  if (!is.null(iso3c) && is.null(contact_matrix_set)) {
    contact_matrix_set <- squire::get_mixing_matrix(iso3c = iso3c)
  }

  country <- get_country(iso3c)

  # Get squire parameters
  pars <- squire::parameters_explicit_SEEIR(
        population = population,
        country = country,
        contact_matrix_set = contact_matrix_set,
        dt = 1, # dt should always be 1 as individual is always a discrete time
        time_period = time_period,
        ...)

  # create list of asymptomatic parameters
  asymp_list <- get_asymptomatic()

  # make final list
  pars <- append(pars, asymp_list)
  pars <- append(pars, list(
    time_period = time_period,
    max_age = max_age
  ))

  # remove all non integer and double class objects as not compatible
  pars <- remove_non_numerics(pars)

  # extend all time varying parameters to fill time steps
  pars$beta_set <- interp_input_par(
    c(pars$tt_beta, time_period),
    c(pars$beta_set, tail(pars$beta_set, 1))
  )

  pars$country <- country

  return(pars)

}


#' @title Append vaccine parameters from nimue model
#' @description This calls \code{\link[nimue]{parameters}} from \href{https://mrc-ide.github.io/nimue/index.html}{nimue}
#' and appends the following parameters to the list:
#'
#'     * \code{gamma_vaccine_delay}: mean duration of period from vaccination to vaccine protection.
#'     * \code{gamma_V}: mean duration of vaccine-derived immunity (days)
#'     * \code{rel_infectiousness_vaccinated}: Relative infectiousness per age category of vaccinated individuals relative to unvaccinated individuals. Default = rep(1, 17), which is no impact of vaccination on onwards transmissions
#'     * \code{rel_infectiousness}: Relative infectiousness per age category relative to maximum infectiousness category. Default = rep(1, 17)
#'     * \code{vaccine_efficacy_infection}: Efficacy of vaccine against infection. This parameter must either be length 1 numeric (a single efficacy for all age groups) or length 17 numeric vector (an efficacy for each age group). An efficacy of 1 will reduce FOI by 100 percent, an efficacy of 0.2 will reduce FOI by 20 percent etc.
#'     * \code{prob_hosp}: Efficacy of vaccine against severe (requiring hospitilisation) disease (by age). This parameter must either be length 1 numeric (a single efficacy for all age groups) or length 17 numeric vector (an efficacy for each age group). An efficacy of 1 will reduce the probability of hospitalisation by 100 percent, an efficacy of 0.2 will reduce the probability of hospitalisation by 20 percent etc.
#'     * \code{vaccine_coverage_mat}: Vaccine coverage targets by age (columns) and priority (row)
#'     * \code{N_prioritisation_steps}: number of vaccine prioritization levels
#'     * \code{current_prioritisation_step}: current prioritization step we are on
#'     * \code{vaccine_set}: vaccines available each day of simulation
#' @param parameters list from [get_parameters]
#' @param ... Other parameters for \code{nimue:::}\code{\link[nimue]{parameters}}, if specifying these parameters, be aware that \code{tt_vaccine_efficacy_infection,tt_vaccine_efficacy_disease} should not be set because safir does not currently model time varying versions of these parameters
#' @export
append_vaccine_nimue <- function(parameters, ...) {

  nimue_pars <- call_nimue_pars(country = parameters$country, ...)

  stopifnot(all(dim(nimue_pars$vaccine_efficacy_infection) == c(1,17,6)))
  stopifnot(all(dim(nimue_pars$prob_hosp) == c(1,17,6)))

  parameters$gamma_vaccine_delay <- 1 / (nimue_pars$gamma_vaccine[2] / 2)
  parameters$gamma_V <- 1 / (nimue_pars$gamma_vaccine[4] / 2)
  parameters$rel_infectiousness_vaccinated <- nimue_pars$rel_infectiousness_vaccinated
  parameters$rel_infectiousness <- nimue_pars$rel_infectiousness
  parameters$vaccine_coverage_mat <- nimue_pars$vaccine_coverage_mat
  parameters$N_prioritisation_steps <- nimue_pars$N_prioritisation_steps
  parameters$current_prioritisation_step <- 1L

  parameters$vaccine_set <- interp_input_par(
    x = c(nimue_pars$tt_vaccine, parameters$time_period),
    y = c(nimue_pars$max_vaccine, tail(nimue_pars$max_vaccine, 1))
  )

  parameters$prob_hosp <- interp_input_par(
    x = c(nimue_pars$tt_vaccine_efficacy_disease, parameters$time_period),
    y = c(nimue_pars$prob_hosp, tail(nimue_pars$prob_hosp, 1))
  )

  # parameters$vaccine_efficacy_infection <- nimue_pars$vaccine_efficacy_infection
  # parameters$prob_hosp <- nimue_pars$prob_hosp

  return(parameters)
}


#' @noRd
call_nimue_pars <- function(
    country,
    # durations
    dur_E  = nimue:::durs$dur_E,
    dur_IMild = nimue:::durs$dur_IMild,
    dur_ICase = nimue:::durs$dur_ICase,
    # hospital durations
    dur_get_ox_survive = nimue:::durs$dur_get_ox_survive,
    dur_get_ox_die = nimue:::durs$dur_get_ox_die,
    dur_not_get_ox_survive = nimue:::durs$dur_not_get_ox_survive,
    dur_not_get_ox_die = nimue:::durs$dur_not_get_ox_die,
    dur_get_mv_survive = nimue:::durs$dur_get_mv_survive,
    dur_get_mv_die = nimue:::durs$dur_get_mv_die,
    dur_not_get_mv_survive = nimue:::durs$dur_not_get_mv_survive,
    dur_not_get_mv_die = nimue:::durs$dur_not_get_mv_die,
    dur_rec = nimue:::durs$dur_rec,
    # vaccine
    dur_R = nimue:::vaccine_pars$dur_R,
    dur_V = nimue:::vaccine_pars$dur_V,
    vaccine_efficacy_infection = nimue:::vaccine_pars$vaccine_efficacy_infection,
    tt_vaccine_efficacy_infection = nimue:::vaccine_pars$tt_vaccine_efficacy_infection,
    vaccine_efficacy_disease = nimue:::vaccine_pars$vaccine_efficacy_disease,
    tt_vaccine_efficacy_disease = nimue:::vaccine_pars$tt_vaccine_efficacy_disease,
    max_vaccine = nimue:::vaccine_pars$max_vaccine,
    tt_vaccine = nimue:::vaccine_pars$tt_vaccine,
    dur_vaccine_delay = nimue:::vaccine_pars$dur_vaccine_delay,
    vaccine_coverage_mat = nimue:::vaccine_pars$vaccine_coverage_mat,
    # health system capacity
    hosp_bed_capacity = NULL,
    ICU_bed_capacity = NULL,
    tt_hosp_beds = 0,
    tt_ICU_beds = 0,
    # misc
    seeding_cases = 20,
    seeding_age_order = NULL
) {
  nimue:::parameters(
    country = country,
    # durations
    dur_E  = dur_E,
    dur_IMild = dur_IMild,
    dur_ICase = dur_ICase,
    # hospital durations
    dur_get_ox_survive = dur_get_ox_survive,
    dur_get_ox_die = dur_get_ox_die,
    dur_not_get_ox_survive = dur_not_get_ox_survive,
    dur_not_get_ox_die = dur_not_get_ox_die,
    dur_get_mv_survive = dur_get_mv_survive,
    dur_get_mv_die = dur_get_mv_die,
    dur_not_get_mv_survive = dur_not_get_mv_survive,
    dur_not_get_mv_die = dur_not_get_mv_die,
    dur_rec = dur_rec,
    # vaccine
    dur_R = dur_R,
    dur_V = dur_V,
    vaccine_efficacy_infection = vaccine_efficacy_infection,
    tt_vaccine_efficacy_infection = tt_vaccine_efficacy_infection,
    vaccine_efficacy_disease = vaccine_efficacy_disease,
    tt_vaccine_efficacy_disease = tt_vaccine_efficacy_disease,
    max_vaccine = max_vaccine,
    tt_vaccine = tt_vaccine,
    dur_vaccine_delay = dur_vaccine_delay,
    vaccine_coverage_mat = vaccine_coverage_mat,
    # health system capacity
    hosp_bed_capacity = hosp_bed_capacity,
    ICU_bed_capacity = ICU_bed_capacity,
    tt_hosp_beds = tt_hosp_beds,
    tt_ICU_beds = tt_ICU_beds,
    # misc
    seeding_cases = seeding_cases,
    seeding_age_order = seeding_age_order
  )
}

#' @title Get population from SQUIRE model
#' @description rounds population sizes to discrete numbers
#'
#' @param iso3c three letter code for your country of interest
#' @export
get_population <- function(iso3c) {
  squire::get_population(iso3c = iso3c, simple_SEIR = FALSE)
}


#' @noRd
get_country <- function(iso3c) {
  squire::population[squire::population$iso3c == iso3c, "country"][[1]]
}


#' Function to add asymptomatic information
#'
#' @return list of asymptomatic parameters
get_asymptomatic <- function() {

  # Temporary values just for testing purposes. Literature derived values
  # will be supplied later
  prob_asymp <- c(0.3, 0.3, rep(0.2, 15))
  IAsymp_0 <- c(rep(0L, 17))
  dur_IAsymp <- squire::default_durations()$dur_IMild

  list(
    prob_asymp = prob_asymp,
    IAsymp_0 = IAsymp_0,
    dur_IAsymp = dur_IAsymp
  )

}


#' @noRd
interp_vector_constant <- function(x, y, by = 1, end = max(x)) {

  # create the time series
  xout <- seq(min(x), end, 1)

  # constant interpolation
  yout <- stats::approx(x,
                        y,
                        xout = xout,
                        method = "constant",
                        yright = tail(y, 1),
                        ties = "ordered")

  return(yout$y)

}


#' @noRd
interp_matrix_list_constant <- function(x, y, by = 1, end = max(x) + 1) {

  # create the time series
  x_all <- c(x, end)
  xout <- seq(min(x), end, 1)

  # build array from matrix list
  rep_mats <- mapply(rep, y, diff(x_all))
  yout <- array(unlist(rep_mats), dim = c(dim(y[[1]]), length(xout)))

  return(yout)

}

# Browse[2]> tt_vaccine_efficacy_disease <<- nimue_pars$tt_vaccine_efficacy_disease
# Browse[2]> prob_hosp <<- nimue_pars$prob_hosp
# Browse[2]> tt_vaccine_efficacy_infection <<- nimue_pars$tt_vaccine_efficacy_infection
# Browse[2]> vaccine_efficacy_infection <<- nimue_pars$vaccine_efficacy_infection


#' @noRd
interp_nimue_array <- function(x, y, by = 1, end = max(x) + 1) {
  stopifnot(inherits(y, "array"))
  stopifnot(length(dim(y)) == 3)
  stopifnot(dim(y)[2] == 17)
  stopifnot(dim(y)[3] == 6)
  stopifnot(dim(y)[1] >= 1)

  # create the time series
  if (length(x) > 1) {
    x_all <- c(x[-length(x)], end)
  } else {
    x_all <- c(x, end)
  }

  xout <- seq(1, end, 1)
  yout <- array(data = NA,dim = c(length(xout),dim(y)[-1]))

  xdiffs <- diff(x_all)
  for (i in 1:6) {
    for (j in seq_along(xdiffs)) {
      yout[(x_all+1)[j]:x_all[j+1], , i] <- do.call(rbind,replicate(n = xdiffs[j],expr = y[j,,i],simplify = FALSE))
    }
  }

}


#' Interpolate input parameters
#'
#' @details Constant interpolation of time changing parameters
#' @param x Time points for \code{y} changing
#' @param y Object to be interpolated
#' @param by Time steps to interpolate at. Default = 1
#' @param end End time point for interpoation. Default = \code{max(x)}
#' @importFrom stats rgamma
interp_input_par <- function(x, y, by = 1, end = max(x)) {

  # check lengths
  assert_length(x, length(y))

  # match input class
  if (inherits(y, "numeric")) {

    ret <- interp_vector_constant(x, y, by = by, end = end)

  } else if (inherits(y, "list") && inherits(y[[1]], "matrix")) {

    ret <- interp_matrix_list_constant(x, y, by = by, end = end)

  } else {

    stop("No interpolation for object y with class ", class(y))

  }

  return(ret)
}
