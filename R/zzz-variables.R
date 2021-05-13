#' @title Create variables
#' @description Create all individual variables for humans
#'
#' @param pop population list
#' @param parameters model parameters
#'
#' @return named list of individual::Variable
#' @noRd
create_variables_zzz <- function(pop, parameters) {

  c(
    create_age_variables(pop, parameters),
    states = create_state_variables(parameters)
  )
}


vaccination_variables_zzz <- function() {

    dose_type <- individual::IntegerVariable$new(discrete_age)

}