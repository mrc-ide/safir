# --------------------------------------------------
#   Translate Ab titre into efficacy
#   and track Ab titre over time
#   Sean L. Wu (slwood89@gmail.com)
#   July 2021
# --------------------------------------------------

#' #' @param ab_titre a \code{\link[individual]{DoubleVariable}} giving the antibody titre for the population
#' #' @param dose_time a list of \code{\link[individual]{IntegerVariable}} objects giving the time of doses for the population
#' #' @param dose_num a \code{\link[individual]{IntegerVariable}} giving the dose number of the population
#' #'
#' vaccine_ab_titre <- function(timestep, dt, ab_titre, dose_time, dose_num) {
#'
#'   # those who have gotten at least one dose
#'   vaccinated <- dose_num$get_index_of(set = 0)
#'   vaccinated$not()
#'
#'   # calculate it....should run each day
#'
#'   # schedule an update
#'
#' }


#' @param ab_titre a \code{\link[individual]{DoubleVariable}} giving the antibody titre for the population
#' @param dose_time a list of \code{\link[individual]{IntegerVariable}} objects giving the time of doses for the population
#' @param dose_num a \code{\link[individual]{IntegerVariable}} giving the dose number of the population
#'
vaccine_ab_titre <- function(parameters, variables, events, dt) {

  dose_num <- variables$dose_num
  dose_time <- variables$dose_time
  ab_titre <- variables$ab_titre

  return(
    function(timestep) {

      # those who have gotten at least one dose
      vaccinated <- dose_num$get_index_of(set = 0)
      vaccinated$not()

      # calculate it....should run each day

      # current Ab titre
      current_ab_values <- ab_titre$get_values(index = vaccinated)

      # apply the discrete difference operator
      # new_ab_values <- blah
      new_ab_values <- pmax(x - (x * 0.01 * dt), 0) # fake dynamics

      # schedule an update
      ab_titre$queue_update(values = new_ab_values, index = vaccinated)

    }
  )

}


#' @param who a \code{\link[individual]{Bitset}} telling who we should calculate efficacy against infection for
vaccine_efficacy_infection <- function(ab_titre, who) {

}


#' @param who a \code{\link[individual]{Bitset}} telling who we should calculate efficacy against severe infection for
vaccine_efficacy_severe <- function(ab_titre, who) {

}
