# --------------------------------------------------
#   additional rendering for safir
#   Sean L. Wu (slwood89@gmail.com)
#   June 2021
# --------------------------------------------------

#' @title Render integer variables
#' @description Renders the number of individuals in each integer count. Should only be used
#' for integer variables with bounded ranges
#' @param renderer a \code{\link[individual]{Render}} object
#' @param variable a \code{\link[individual]{IntegerVariable}} object
#' @param margin the set of values which may be taken for the integer variable
#' @export
integer_count_render_process <- function(renderer, variable, margin) {
  stopifnot(inherits(variable, "IntegerVariable"))
  stopifnot(inherits(renderer, "Render"))
  function(t) {
    for (c in margin) {
      renderer$render(paste0(c, '_count'), variable$get_size_of(c), t)
    }
  }
}

#' @title Render integer variables every day
#' @description Renders the number of individuals in each integer count. Should only be used
#' for integer variables with bounded ranges. This only renders output on timesteps
#' that correspond to a day. This saves memory and potentially time for
#' simulations with small timesteps.
#' @param renderer a \code{\link[individual]{Render}} object
#' @param variable a \code{\link[individual]{IntegerVariable}} object
#' @param margin the set of values which may be taken for the integer variable
#' @param dt size of time step
#' @export
integer_count_render_process_daily <- function(renderer, variable, margin, dt) {
  stopifnot(inherits(variable, "IntegerVariable"))
  stopifnot(inherits(renderer, "Render"))
  function(t) {
    if ((t * dt) %% 1 == 0) {
      for (c in margin) {
        renderer$render(paste0(c, '_count'), variable$get_size_of(c), as.integer(t * dt))
      }
    }
  }
}

#' @title Render doses by age every day
#' @description Renders the number of individuals with each dose in each age bin
#' (cross tabulated doses and ages).
#' This only renders output on timesteps that correspond to a day.
#' See examples for how to quickly summarize the output.
#' @param renderer a \code{\link[individual]{Render}} object
#' @param age a \code{\link[individual]{IntegerVariable}} object
#' @param dose a \code{\link[individual]{IntegerVariable}} object
#' @param parameters model parameters
#' @param dt size of time step
#' @examples
#' # if the rendered is called dose_age_renderer
#' tmp <- as.data.table(dose_age_renderer$to_dataframe())
#' tmp <- melt(tmp,id.vars="timestep")
#' tmp1 <- tmp[, tstrsplit(variable, "_", keep = c(2, 4))]
#' tmp[ , variable := NULL]
#' tmp <- cbind(tmp, tmp1)
#' setnames(x = tmp,old = c("V1", "V2"),new = c("dose","age"))
#'
#' dose_dt <- tmp[, .(value = sum(value)), by = .(dose, timestep)]
#' @export
dose_age_render_process_daily <- function(renderer, age, dose, parameters, dt) {
  stopifnot(inherits(age, "IntegerVariable"))
  stopifnot(inherits(dose, "IntegerVariable"))
  stopifnot(inherits(renderer, "Render"))
  num_dose <- parameters$N_phase
  num_age <- parameters$N_age
  function(t) {
    if ((t * dt) %% 1 == 0) {
      day <- as.integer(t * dt)
      dose_age_tab <- cross_tab_doses_age(doses = dose$.variable, age = age$.variable, num_doses = num_dose, num_ages = num_age)
      for (d in 1:nrow(dose_age_tab)) {
        for (a in 1:ncol(dose_age_tab)) {
          renderer$render(paste0("dose_", d-1, "_age_", a, "_count"), dose_age_tab[d, a], day)
        }
      }
    }
  }
}

#' @title Render categories every day
#' @description Renders the number of individuals in each category, but only
#' on timesteps that correspond to a day. This saves memory and potentially time
#' for simulations with small timesteps.
#' @param renderer a \code{\link[individual]{Render}} object
#' @param variable a \code{\link[individual]{CategoricalVariable}} object
#' @param categories a character vector of categories to render
#' @param dt size of time step
#' @return a function which can be passed as a process to \code{\link{simulation_loop}}
#' @export
categorical_count_renderer_process_daily <- function(renderer, variable, categories, dt) {
  stopifnot(inherits(variable, "CategoricalVariable"))
  stopifnot(inherits(renderer, "Render"))
  function(t) {
    if ((t * dt) %% 1 == 0) {
      for (c in categories) {
        renderer$render(paste0(c, '_count'), variable$get_size_of(c), as.integer(t * dt))
      }
    }
  }
}
