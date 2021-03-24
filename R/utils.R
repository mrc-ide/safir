`%||%` <- function(a, b) { # nolint
  if (is.null(a)) b else a
}

#' Erlang waiting time distribution
#'
#' @details Random draws from erlang distribution
#' @param size Number of draws
#' @param mu Mean duration
#' @importFrom stats rgamma
r_erlang <- function(size, mu) { rgamma(size, shape = 2, rate = 2 / mu) }

#' Exponential waiting time distribution
#'
#' @details Random draws from exponential distribution
#' @param size Number of draws
#' @param mu Mean duration
#' @importFrom stats rexp
r_exp <- function(size, mu) { rexp(size, rate = 1 / mu) }

#' @noRd
remove_non_numerics <- function(l) {
  clean <- list()
  for (key in names(l)) {
    if (storage.mode(l[[key]]) %in% c("integer", "double")) {
      clean[[key]] <- l[[key]]
    }
  }
  clean
}

#' @title Bernoulli sample
#' @param p vector of probabilites for each individual
#' @noRd
#' @importFrom stats runif
bernoulli_multi_p <- function(p) runif(length(p), 0, 1) < p

#' @noRd
vcapply <- function(X, FUN, ...) vapply(X, FUN, ..., character(1))

#'@noRd
rep_list <- function(x, times) {
  setNames(lapply(x, function(x) {rep(list(x), times)}), names(x))
}
