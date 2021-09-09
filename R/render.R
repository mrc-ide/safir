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
