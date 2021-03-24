#' @title Define the human model
#' @description Declares the human individual and assigns the
#' relevant states and variables
#'
#' @noRd
#' @param states available states to assign
#' @param variables available variables to assign
#' @param events available events to assign
create_human <- function(states, variables, events) {
  individual::Individual$new(
    "human",
    states = states,
    variables = variables,
    events = events
  )
}
