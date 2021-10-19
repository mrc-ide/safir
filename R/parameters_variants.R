# --------------------------------------------------------------------------------
#   parameters for variants
#   Sean L. Wu (slwood89@gmail.com)
#   August 2021
# --------------------------------------------------------------------------------


#' @title Get variant of concern transmission parameters
#' @param safir_parameters either a single object from \code{\link{get_parameters}}
#' or a list of objects from \code{\link{get_parameters}}
#' @param voc_types character vector of variant names (must include "wt" and "none")
#' @param voc_trajectory
#' @param vaccine_parameters either a single object from \code{\link{get_vaccine_ab_titre_parameters}}
#' or a list of objects from \code{\link{get_vaccine_ab_titre_parameters}} giving efficacy parameters
#' for each variant
#' @export
get_voc_parameters <- function(
  safir_parameters,
  voc_types = c("alpha", "beta", "delta", "wt"),
  voc_trajectory = NULL,
  vaccine_parameters
) {

  stopifnot(is.character(voc_types))
  stopifnot("wt" %in% voc_types)

  parameters <- list()
  parameters$voc_types <- voc_types
  parameters$voc_num <- length(voc_types)


  if (!is.null(attr(safir_parameters, "type"))) {

    stopifnot(attr(safir_parameters, "type") == "safir_squire")
    tmax <- safir_parameters$time_period

  } else if (length(safir_parameters) == parameters$voc_num) {

    stopifnot( all(vapply(X = safir_parameters, FUN = function(l){attr(l, "type")}, FUN.VALUE = character(1)) == "safir_squire") )
    tmax <- safir_parameters[[1]]$time_period

  } else {
    stop("invalid object passed for safir_parameters")
  }

  stopifnot(is.finite(tmax))
  stopifnot(tmax > 1)

  if (is.null(voc_trajectory)) {
    voc_trajectory <- matrix(data = 0, nrow = tmax,ncol = length(voc_types),dimnames = list(NULL, voc_types))
    voc_trajectory[, "wt"] <- 1
  }

  stopifnot(all(rowSums(voc_trajectory) == 1))
  stopifnot(ncol(voc_trajectory) == parameters$voc_num)
  stopifnot(nrow(voc_trajectory) == tmax)

  parameters$voc_trajectory <- voc_trajectory

  if (!is.null(attr(vaccine_parameters, "type"))) {

    stopifnot(attr(vaccine_parameters, "type") == "ab_titre")

    parameters$dr_vec <- vaccine_parameters$dr_vec
    parameters$mu_ab <- vaccine_parameters$mu_ab
    parameters$std10 <- vaccine_parameters$std10
    parameters$correlated <- vaccine_parameters$correlated

    parameters$k <- rep(vaccine_parameters$k, parameters$voc_num)
    parameters$ab_50 <- rep(vaccine_parameters$ab_50, parameters$voc_num)
    parameters$ab_50_severe <- rep(vaccine_parameters$ab_50_severe, parameters$voc_num)

  } else if (length(vaccine_parameters) == parameters$voc_num) {

    stopifnot( all(vapply(X = vaccine_parameters, FUN = function(l){attr(l, "type")}, FUN.VALUE = character(1)) == "ab_titre") )

    parameters$dr_vec <- vaccine_parameters[[1]]$dr_vec
    parameters$mu_ab <- vaccine_parameters[[1]]$mu_ab
    parameters$std10 <- vaccine_parameters[[1]]$std10
    parameters$correlated <- vaccine_parameters[[1]]$correlated

    parameters$ab_50 <- vapply(X = vaccine_parameters, FUN = function(l){l$ab_50}, FUN.VALUE = numeric(1))
    parameters$ab_50_severe <- vapply(X = vaccine_parameters, FUN = function(l){l$ab_50_severe}, FUN.VALUE = numeric(1))
    parameters$k <- vapply(X = vaccine_parameters, FUN = function(l){l$k}, FUN.VALUE = numeric(1))

  } else {
    stop("invalid object passed for vaccine_parameters")
  }

  attr(parameters, "type") <- "safir_vaccine_voc"
  return(parameters)
}
