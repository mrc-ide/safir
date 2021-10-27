#' @title Create natural immunity variables
#' @param variables a list
#' @param parameters list of model parameters
#' @importFrom individual IntegerVariable
#' @return named list of individual::Variable
#' @export
create_natural_immunity_variables <- function(variables, parameters) {

  n <- sum(parameters$population)

  variables$inf_num <- IntegerVariable$new(initial_values = rep(0L, n))
  variables$inf_time <- IntegerVariable$new(initial_values = rep(-1L, n))

  return(variables)
}
