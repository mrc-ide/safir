# --------------------------------------------------------------------------------
#   parameters for variants
#   Sean L. Wu (slwood89@gmail.com)
#   August 2021
# --------------------------------------------------------------------------------


#' @title Get variant of concern parameters
#' @param tmax the maximum day of simulation (not the number of time steps)
#' @param voc_types character vector of variant names (must include "wt" and "none")
#' @export
get_voc_parameters <- function(
  tmax,
  voc_types = c("alpha", "beta", "delta", "wt", "none"),
  voc_trajectory = NULL
) {

  stopifnot(is.finite(tmax))
  stopifnot(tmax > 1)

  stopifnot(is.character(voc_types))
  stopifnot(c("wt","none") %in% voc_types)

  if (is.null(voc_trajectory)) {
    vocs <- voc_types[voc_types != "none"]
    voc_trajectory <- matrix(data = 0, nrow = tmax,ncol = length(vocs),dimnames = list(NULL, vocs))
    voc_trajectory[, "wt"] <- 1
  }

  stopifnot(all(rowSums(voc_trajectory) == 1))
  stopifnot(ncol(voc_trajectory) == (length(voc_types) - 1))
  stopifnot(nrow(voc_trajectory) == tmax)

  parameters <- list()
  parameters$voc_types <- voc_types
  parameters$voc_trajectory <- voc_trajectory

  return(parameters)
}
