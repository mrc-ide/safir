# --------------------------------------------------
#   model parameters for squire transmission model
#   Sean L. Wu (slwood89@gmail.com)
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
