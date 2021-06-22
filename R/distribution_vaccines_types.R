# --------------------------------------------------
#   functions to work with scheduling/allocating vaccines
#   for the general vaccine model but with types
#   Sean L. Wu (slwood89@gmail.com)
#   June 2021
# --------------------------------------------------


#' @title Get proportion of an age group that is vaccinated (multi-dose, with types)
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
  vaccinated_bset <- variables$dose_time[[dose]]$get_index_of(set = -1) # haven't gotten this dose
  vaccinated_bset <- vaccinated_bset$not() # complement = vaccinated people
  vaccinated_bset$and(age_bset)$and(type_bset)
  return( vaccinated_bset$size() / N )
}
