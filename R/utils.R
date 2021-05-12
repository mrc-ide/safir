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

#' Make geometric waiting time distribution
#'
#' @details Random draws from a geometric distribution
#' which approximates a continuous time exponential distribution.
#' See [pomp documentation](https://kingaa.github.io/pomp/vignettes/C_API.html#exponentialgeometric-rate-conversion)
#' for details on the approximation.
#' @param mu Mean duration
#' @param dt Size of time step
#' @importFrom stats rgeom pexp
make_rexp <- function(mu, dt) {
  assert_pos(mu, zero_allowed = FALSE)
  assert_pos(dt, zero_allowed = FALSE)
  R <- 1 / mu
  r <- log(1 + (R*dt)) / dt
  function(n) {
    rgeom(n = n, prob = pexp(q = r * dt)) + 1L
  }
}

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
