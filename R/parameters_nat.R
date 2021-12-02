#' @title Make vector of variant fold reduction in natural antibody titre values
#' @description Generate a vector for each time step where no VFR effect corresponds
#' to values in the vector of `1`, and the appropriate VFR adjustment otherwise.
#' @param parameters a list of model parameters
#' @param dt time step size
#' @param vfr a positive value giving the fold reduction
#' @param vfr_time_1 time to begin reduction
#' @param vfr_time_2 time to end reduction
#' @export
variant_fold_reduction_vector <- function(parameters, dt, vfr, vfr_time_1, vfr_time_2) {

  stopifnot(is.finite(dt))

  stopifnot(is.finite(vfr))
  stopifnot(is.finite(vfr_time_1))
  stopifnot(is.finite(vfr_time_2))

  stopifnot(vfr > 0)
  stopifnot(vfr_time_1 > 0)
  stopifnot(vfr_time_2 > 0)
  stopifnot(vfr_time_2 > vfr_time_1)
  stopifnot(vfr_time_2 < parameters$time_period)

  timesteps <- parameters$time_period / dt
  vfr_time_1_dt <- vfr_time_1 / dt
  vfr_time_2_dt <- vfr_time_2 / dt

  vfr_vec <- c(rep(1, times = vfr_time_1_dt - 1), seq(from = 1, to = vfr, length.out = (vfr_time_2_dt - vfr_time_1_dt + 1)), rep(vfr, times = timesteps - vfr_time_2_dt))

  return(vfr_vec)
}


#' @title Attach parameters dealing with immune (NAT) dynamics to parameter list
#' @param parameters a list of model parameters
#' @param vfr a vector from [safir::variant_fold_reduction_vector]
#' @param mu_ab_infection a vector of values associated with NAT for each infection bout;
#' please read [safir::attach_event_listeners_natural_immunity] carefully depending on
#' if you are using the additive or overwriting NAT model, because this parameter's
#' interpretation will change
#' @param std10_infection a standard deviation for log-normal draw of NAT values
#' @export
make_immune_parameters <- function(parameters, vfr, mu_ab_infection = NULL, std10_infection = NULL) {

  if (!is.null(std10_infection)) {
    stopifnot(length(std10_infection) == 1L)
    stopifnot(is.finite(std10_infection))
    stopifnot(std10_infection > 0)
  }

  if (!is.null(mu_ab_infection)) {
    stopifnot(length(mu_ab_infection) > 0)
    stopifnot(is.finite(mu_ab_infection))
    stopifnot(mu_ab_infection > 0)
  }

  stopifnot(length(vfr) >= parameters$time_period)

  if (is.null(mu_ab_infection)) {
    stopifnot(!is.null(parameters$mu_ab))
    mu_ab_infection <- parameters$mu_ab
  }

  if (is.null(std10_infection)) {
    stopifnot(!is.null(parameters$std10))
    std10_infection <- parameters$std10
  }

  parameters$vfr <- vfr
  parameters$mu_ab_infection <- mu_ab_infection
  parameters$std10_infection <- std10_infection

  return(parameters)
}
