# --------------------------------------------------
#   additional rendering for safir
#   Sean L. Wu (slwood89@gmail.com)
#   June 2021
# --------------------------------------------------

# should prob move to individual eventually


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
