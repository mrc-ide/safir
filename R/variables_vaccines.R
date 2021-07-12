# --------------------------------------------------
#   Create and query vaccination variables
#   Sean L. Wu (slwood89@gmail.com)
#   June 2021
# --------------------------------------------------


# --------------------------------------------------
#   nimue version
# --------------------------------------------------

#' @title Create vaccination variables (nimue vaccine model)
#' @description Create all individual variables for humans
#'
#' @param variables a list
#' @param pop population list
#'
#' @return named list of individual::Variable
#' @export
create_vaccine_variables_nimue <- function(variables, pop) {

  n <- sum(pop$n)

  variables$vaccine_states <- individual::IntegerVariable$new(initial_values = rep(1,n))
  variables$eligible <- individual::Bitset$new(size = n)
  variables$vaccinated <- individual::Bitset$new(size = n)
  variables$empty <- individual::Bitset$new(size = n)

  return(variables)
}

#' @title Get proportion of an age group that is vaccinated (nimue vaccine model)
#' @description Get proportion of an age group that has received a particular vaccine dose
#' by this timestep.
#' @param variables a list
#' @param age an age group (can be a vector of multiple age groups)
#'
#' @export
get_proportion_vaccinated_nimue <- function(variables, age) {
  age_bset <- variables$discrete_age$get_index_of(age)
  N <- age_bset$size()
  vaccinated_bset <- variables$vaccinated$copy()
  vaccinated_bset$and(age_bset)
  return( vaccinated_bset$size() / N )
}


# --------------------------------------------------
#   multiple doses, no types
# --------------------------------------------------

#' @title Create vaccination variables (multi-dose, no types)
#' @description Create all individual variables for humans
#'
#' @param variables a list
#' @param pop population list
#' @param max_dose maximum number of possible doses
#'
#' @return named list of individual::Variable
#' @export
create_vaccine_variables <- function(variables, pop, max_dose = 2) {

  n <- sum(pop)

  variables$ab_titre <- individual::DoubleVariable$new(initial_values = rep(0,n))
  variables$dose_num <- individual::IntegerVariable$new(initial_values = rep(0,n))
  # dose time stores the time step when it happened, not necessarily the same as the day (if dt != 1)
  variables$dose_time <- replicate(n = max_dose,expr = individual::IntegerVariable$new(initial_values = rep(-1,n)),simplify = FALSE)
  variables$phase <- new.env(hash = FALSE)
  variables$phase$value <- 1L

  return(variables)
}

#' @title Initialize vaccination variables (multi-dose, no types)
#' @param variables a list from \code{\link{create_vaccine_variables}}
#' @param dose_time_init list of dose times (must all be integer vectors of the same length, -1 is value to indicate that individual has not
#' yet received this dose)
#' @param dose_num_init vector of initial doses (must be self-consistent with \code{dose_time_init})
#' @export
initialize_vaccine_variables <- function(variables, dose_time_init, dose_num_init) {
  stopifnot(inherits(dose_time_init, "list"))
  stopifnot(length(dose_time_init) > 0)
  stopifnot(length(dose_time_init[[1]]) == length(dose_num_init))
  # if you haven't gotten the first dose, you are at dose 0
  stopifnot(all(which(dose_time_init[[1]] < 0) == which(dose_num_init == 0)))
  for (d in seq_along(dose_time_init)) {
    variables$dose_time[[d]]$queue_update(values = dose_time_init[[d]])
    variables$dose_time[[d]]$.update()
  }
  variables$dose_num$queue_update(value = dose_num_init)
  variables$dose_num$.update()
}

#' @title Update vaccine variables
#' @description This should be called from the simulation loop. It does not
#' update disease state.
#' @param variables a list
#' @export
update_vaccine_variables <- function(variables) {

  for (dose in variables$dose_time) {
    dose$.update()
  }
  variables$dose_num$.update()

}
