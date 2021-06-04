# --------------------------------------------------
#   utility functions
#   May 2021
#   1. `%||%`
#   2. make_rerlang
#   3. make_rexp
#   4. make_rexp_simple
#   5. remove_non_numerics
#   6. vcapply
#   7. rep_list
# --------------------------------------------------


`%||%` <- function(a, b) { # nolint
  if (is.null(a)) b else a
}


#' @title Make Erlang waiting time distribution
#'
#' @description Random draws from Erlang distribution
#' scaled by time step size.
#' @param mu Mean duration
#' @param dt Size of time step
#' @param shape Shape parameter of Erlang distribution
#' @param shift number of time steps to add to sampled value
#' @importFrom stats rgamma
#' @export
make_rerlang <- function(mu, dt, shape = 2, shift = 0L) {
  assert_pos(mu, zero_allowed = FALSE)
  assert_pos(dt, zero_allowed = FALSE)
  assert_pos(shape, zero_allowed = FALSE)
  r <- shape / mu
  function(n) {
    # x <- floor(rgamma(n = n, shape = shape, rate = r) / dt) + shift # possibility of 0 delay
    # if (any(is.na(x))) {
    #   browser()
    # }
    # return(x)
    floor(rgamma(n = n, shape = shape, rate = r) / dt) + shift # possibility of 0 delay
  }
}


#' @title Make discretized exponential waiting time distribution
#'
#' @description Random draws from a geometric distribution
#' which approximates a continuous time exponential distribution.
#' See [pomp documentation](https://kingaa.github.io/pomp/vignettes/C_API.html#exponentialgeometric-rate-conversion)
#' for details on the approximation.
#' @param mu Mean duration
#' @param dt Size of time step
#' @param shift number of time steps to add to sampled value
#' @importFrom stats rgeom pexp
make_rexp <- function(mu, dt, shift = 0L) {
  assert_pos(mu, zero_allowed = FALSE)
  assert_pos(dt, zero_allowed = FALSE)
  R <- 1 / mu
  r <- log(1 + (R*dt)) / dt
  p <- pexp(q = r * dt)
  function(n) {
    rgeom(n = n, prob = p) + shift
  }
}


#' @title Make discretized exponential waiting time distribution (simple)
#' @description Make simple geometric approximation to
#' continuous time exponential distribution.
#' @param mu Mean duration
#' @param dt Size of time step
#' @param shift number of time steps to add to sampled value
#' @importFrom stats rgeom pexp
make_rexp_simple <- function(mu, dt, shift = 0L) {
  assert_pos(mu, zero_allowed = FALSE)
  assert_pos(dt, zero_allowed = FALSE)
  r <- 1 / mu
  p <- pexp(q = r * dt)
  function(n) {
    rgeom(n = n, prob = p) + shift
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


#' @noRd
vcapply <- function(X, FUN, ...) vapply(X, FUN, ..., character(1))


#'@noRd
rep_list <- function(x, times) {
  setNames(lapply(x, function(x) {rep(list(x), times)}), names(x))
}
