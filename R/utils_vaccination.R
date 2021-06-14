# --------------------------------------------------
#   functions to work with scheduling/allocating vaccines
#   for the general vaccine model
#   Sean L. Wu (slwood89@gmail.com)
#   June 2021
# --------------------------------------------------

#' @title Get proportion of an age group that is vaccinated
#' @description Get proportion of an age group that has received a particular vaccine dose
#' by this timestep. This is similar to the function \code{\link[nimue]{coverage}} in the nimue package.
#' @param variables a list
#' @param age an age group (can be a vector of Multiple age groups)
#' @param dose which dose to get proportion of age group who has received it
#'
#' @export
get_proportion_vaccinated <- function(variables, age, dose) {
  age_bset <- variables$discrete_age$get_index_of(age)
  N <- age_bset$size()
  vaccinated_bset <- variables$dose_time[[dose]]$get_index_of(set = -1) # not vaccinated people
  vaccinated_bset <- vaccinated_bset$not() # complement = vaccinated people
  vaccinated_bset$and(age_bset)
  return( vaccinated_bset$size() / N )
}

#' @title Get proportion of an age group that is vaccinated, by type
#' @description Get proportion of an age group that has received a particular vaccine type and dose
#' by this timestep. This is similar to the function \code{\link[nimue]{coverage}} in the nimue package.
#' @param variables a list
#' @param age an age group (can be a vector of Multiple age groups)
#' @param type type of vaccine
#' @param dose which dose to get proportion of age group who has received it
#'
#' @export
get_proportion_vaccinated_type <- function(variables, age, type, dose) {
  # stopifnot(type %in% variables$dose_type[[dose]]$get_categories())
  age_bset <- variables$discrete_age$get_index_of(age)
  N <- age_bset$size()
  type_bset <- variables$dose_type[[dose]]$get_index_of(type)
  vaccinated_bset <- variables$dose_time[[dose]]$get_index_of(set = -1) # not vaccinated people
  vaccinated_bset <- vaccinated_bset$not() # complement = vaccinated people
  vaccinated_bset$and(age_bset)$and(type_bset)
  return( vaccinated_bset$size() / N )
}

#' @title Identity those persons eligible for a dose
#' @description Find those individuals who have had the dose preceding \code{dose_number},
#' have not yet received the next one, and are beyond the \code{dose_period}.
#' This is similar to the function \code{\link[nimue]{eligable_for_second}} in the nimue package.
#' @param dose_number which dose? (must be greater than 1)
#' @param dose_period days between \code{dose_number} and the previous dose (please make sure \code{dose_number / dt} produces an integer number of timesteps)
#' @param variables a list
#' @param t the current time step
#' @param dt size of the time step
#' @return an \code{\link[individual]{Bitset}}
#' @export
eligible_for_dose_vaccine <- function(dose_number, dose_period, variables, t, dt) {
  stopifnot(dose_number > 1)
  # who has gotten the previous dose? (with correction for dt < 1)
  had_previous_beyond_threshold <- variables$dose_time[[dose_number - 1]]$get_index_of(a = 0, b = t - as.integer(dose_period/dt))
  # who has not gotten the next one?
  not_had_next <- variables$dose_time[[dose_number]]$get_index_of(set = -1)
  # return people past the threshold and who haven't gotten the next one yet
  return(not_had_next$and(had_previous_beyond_threshold))
}
