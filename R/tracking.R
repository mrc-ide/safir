# --------------------------------------------------
#   tracking output
#   Sean L. Wu (slwood89@gmail.com)
#   September 2021
# --------------------------------------------------

#' @title Create listener to track incidence
#'
#' @param renderer a \code{\link[individual]{Render}} object
#' @export
create_incidence_tracking_listener <- function(renderer) {
  function(timestep, target) {
    renderer$render(name = "incidence", value = target$size(), timestep = timestep)
  }
}

#' @title Attach listener to track incidence
#' @description To use incidence tracking, run the lines in examples
#' anytime after all events have been made and initialized and before
#' simulation begins.
#' @param events a list of events
#' @param renderer a \code{\link[individual]{Render}} object
#' @examples
#' \dontrun{
#' incidence_renderer <- individual::Render$new(timesteps)
#' attach_tracking_listener_incidence(events = events,renderer = incidence_renderer)
#' }
#' @export
attach_tracking_listener_incidence <- function(events, renderer) {
  events$exposure$add_listener(
    create_incidence_tracking_listener(renderer = renderer)
  )
}
