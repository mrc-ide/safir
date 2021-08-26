# --------------------------------------------------
#   Create and query variant variables
#   Sean L. Wu (slwood89@gmail.com)
#   August 2021
# --------------------------------------------------


#' @title Create variant of concern variables
#' @description Please note the variant model assumes the user is applying
#' the vaccine model with antibody dynamics (not nimue or squire versions).
#' @param variables a list
#' @param parameters list of model parameters
#' @importFrom individual CategoricalVariable
#' @return named list of individual::Variable
#' @export
create_voc_variables <- function(variables, parameters) {

  stopifnot(all(c("wt","none") %in% parameters$voc_types))
  n <- sum(parameters$population)

  # this tracks what a person's current infection is
  variables$voc <- CategoricalVariable$new(categories = parameters$voc_types,initial_values = rep("none", n))

  return(variables)
}

#' @title Update variant of concern variables
#' @description This should be called from the simulation loop. It does not
#' update disease state.
#' @param variables a list
#' @export
update_voc_variables <- function(variables) {

  variables$voc$.update()

}
