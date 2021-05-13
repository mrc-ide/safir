# --------------------------------------------------
#   utility functions
#   May 2021
#   1. make_rerlang
#   2. make_rexp
# --------------------------------------------------

#' Make Erlang waiting time distribution
#'
#' @details Random draws from Erlang distribution
#' scaled by time step size.
#' @param mu Mean duration
#' @param dt Size of time step
#' @param shape Shape parameter of Erlang distribution
#' @importFrom stats rnbinom
make_rerlang <- function(mu, dt, shape = 2) {
  assert_pos(mu, zero_allowed = FALSE)
  assert_pos(dt, zero_allowed = FALSE)
  assert_pos(shape, zero_allowed = FALSE)
  r <- shape / mu
  function(n) {
    rgamma(n = n, shape = shape, rate = r) / dt # possibility of 0 delay
  }
}

#' Make discretized exponential waiting time distribution
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
  p <- pexp(q = r * dt)
  function(n) {
    rgeom(n = n, prob = p) + 1L
  }
}
