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
#   full version
# --------------------------------------------------

#' @title Create vaccination variables
#' @description Create all individual variables for humans
#'
#' @param variables a list
#' @param pop population list
#' @param vaxx_types a character vector of possible vaccine types
#' @param max_dose maximum number of possible doses
#'
#' @return named list of individual::Variable
#' @export
create_vaccine_variables <- function(variables, pop, vaxx_types = c("AZ", "MD", "PF", "JJ", "novax"), max_dose = 2) {

  stopifnot("novax" %in% vaxx_types)
  stopifnot(length(novax) >= 2)

  n <- sum(pop$n)

  variables$dose_num <- individual::IntegerVariable$new(initial_values = rep(0,n))
  variables$dose_time <- replicate(n = max_dose,expr = individual::IntegerVariable$new(initial_values = rep(-1,n)),simplify = FALSE)
  variables$dose_type <- replicate(n = max_dose,expr = individual::CategoricalVariable$new(categories = vaxx_types,initial_values = rep("novax",n)),simplify = FALSE)
  variables$phase <- 1L

  return(variables)
}



