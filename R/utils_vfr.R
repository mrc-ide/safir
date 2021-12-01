#' @title Make vector of variant fold reduction in natural antibody titre values
#' @description Generate a vector for each time step where no VFR effect corresponds
#' to values in the vector of `1`, and the appropriate VFR adjustment otherwise.
#' @param parameters a list of model parameters
#' @param variables a list of model variables
#' @param dt time step size
#' @param vfr a list of model parameters
#' @param vfr_time_1 a list of model variables
#' @param vfr_time_2 time step size
#' @export
variant_fold_reduction_vector <- function(parameters, variables, dt, vfr, vfr_time_1, vfr_time_2) {

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
