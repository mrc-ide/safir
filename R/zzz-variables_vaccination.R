# --------------------------------------------------
#   vaccination variables and helper fn
#   May 2021
#   1. create_vaccine_variables
#   2. get_proportion_vaccinated
# --------------------------------------------------


# --------------------------------------------------
#   nimue version
# --------------------------------------------------

#' @title Create vaccination variables for nimue style vaccinations
#' @description Create all individual variables for humans
#'
#' @param variables a list
#' @param pop population list
#'
#' @return named list of individual::Variable
#' @export
create_vaccine_variables_nimue <- function(variables, pop) {

  n <- sum(pop$n)

  variables$vaccine_state <- individual::CategoricalVariable$new(categories = c("v0","v1v2","v3v4","v5"),initial_values = rep("v0",n))
  variables$time <- individual::IntegerVariable$new(initial_values = rep(-1,n))
  variables$vaccinated <- individual::Bitset$new(size = n)
  variables$eligible <- individual::Bitset$new(size = n)
  variables$empty <- individual::Bitset$new(size = n)

  return(variables)
}

#' @title Get proportion of an age group that is vaccinated for nimue style vaccination
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
#   full version
# --------------------------------------------------

#' @title Create vaccination variables
#' @description Create all individual variables for humans
#'
#' @param pop population list
#' @param vaxx_types a character vector of possible vaccine types
#' @param max_dose maximum number of possible doses (2 for now)
#'
#' @return named list of individual::Variable
#' @export
create_vaccine_variables <- function(pop, vaxx_types = c("AZ", "MD", "PF", "JJ", "null"), max_dose = 2) {

  n <- sum(pop$n)
  vaxx <- list()

  vaxx$dose_time <- replicate(n = max_dose,expr = individual::IntegerVariable$new(initial_values = rep(-1,n)),simplify = FALSE)
  vaxx$dose_type <- replicate(n = max_dose,expr = individual::CategoricalVariable$new(categories = vaxx_types,initial_values = rep("null",n)))

  return(vaxx)
}


#' @title Get proportion of an age group that is vaccinated
#' @description Get proportion of an age group that has received a particular vaccine dose
#' by this timestep.
#' @param variables a list
#' @param timestep the current time step
#' @param age an age group (can be a vector of Multiple age groups)
#' @param dose which dose to get proportion of age group who has received it
#'
#' @export
get_proportion_vaccinated <- function(variables, timestep, age, dose) {
  age_bset <- variables$discrete_age$get_index_of(age)
  N <- age_bset$size()
  vaccinated_bset <- variables$dose_time[[dose]]$get_index_of(a = 0,b = timestep)
  vaccinated_bset$and(age_bset)
  return( vaccinated_bset$size() / N )
}

