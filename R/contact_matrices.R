#' @title Contact matrix at current timestep
#' @description Currently only supports one contact matrix
#' @param contact_matrix_set model parameter for the contact matrices
#' @return 2D contact matrix for the current timestep
#' @noRd
get_contact_matrix <- function(contact_matrix_set) {
  contact_matrix_set[1, , ]
}
