# --------------------------------------------------
#   contact matrix
#   May 2021
#   1. get_contact_matrix
# --------------------------------------------------


#' @title Contact matrix at current timestep
#' @description Currently only supports one contact matrix
#' @param parameters object returned from [get_parameters]
#' @return 2D contact matrix for the current timestep
#' @noRd
get_contact_matrix <- function(parameters) {
  parameters$mix_mat_set[1, , ]
}
