# --------------------------------------------------
#   events to handle modeling of natural immunity
#   Sean L. Wu (slwood89@gmail.com)
#   October 2021
# --------------------------------------------------

#' @title Attach event listeners for natural immunity events
#' @param variables a named list of variables
#' @param events a named list of events
#' @param parameters the parameters
#' @param dt size of time step
#' @param additive if `FALSE` the antibody titre is overwritten upon recovery from infection
#' according to `parameters$mu_ab_infection`, if `TRUE` the antibody titre is sampled and added
#' to existing titre. Please choose this option with caution, as it changes the interpretation
#' of the parameter value `parameters$mu_ab_infection`.
#' @export
attach_event_listeners_natural_immunity <- function(variables, events, parameters, dt, additive = FALSE) {

  # checks
  stopifnot(is.logical(additive))

  stopifnot(!is.null(parameters$mu_ab_infection))
  stopifnot(is.finite(parameters$mu_ab_infection))

  if (is.null(parameters$std10_infection)) {
    std10_infection <- parameters$std10
  } else {
    stopifnot(length(parameters$std10_infection) == 1L)
    stopifnot(is.finite(parameters$std10_infection))
    std10_infection <- parameters$std10_infection
  }

  # recovery: handle 1 timestep R->S and update ab titre for immune response
  if (length(events$recovery$.listeners) == 2) {
    events$recovery$.listeners[[2]] <- NULL
  }

  # they go from R to S in 1 time step
  events$recovery$add_listener(
    function(timestep, target) {
      events$immunity_loss$schedule(target = target, delay = rep(1, target$size()))
    }
  )

  # effect of infection on NAT
  if (additive) {
    events$recovery$add_listener(
      function(timestep, target) {

        # update inf_num
        inf <- variables$inf_num$get_values(target) + 1L
        variables$inf_num$queue_update(values = inf, index = target)

        # update last time of infection
        variables$inf_time$queue_update(values = timestep, index = target)

        # get NAT values and convert to linear scale
        current_ab_titre <- variables$ab_titre$get_values(index = target)
        current_ab_titre <- exp(current_ab_titre)

        # draw NAT boost on linear scale
        inf[inf > length(parameters$mu_ab_infection)] <- length(parameters$mu_ab_infection)
        zdose <- 10^rnorm(n = target$size(), mean = log10(parameters$mu_ab_infection[inf]),sd = std10_infection)
        new_ab_titre <- current_ab_titre + zdose

        # back to ln scale, and impose max value constraint
        new_ab_titre <- log(new_ab_titre)
        new_ab_titre <- pmin(new_ab_titre, parameters$max_ab)

        # queue NAT update
        variables$ab_titre$queue_update(values = new_ab_titre, index = target)
      }
    )
  } else {
    events$recovery$add_listener(
      function(timestep, target) {
        # update inf_num
        inf <- variables$inf_num$get_values(target) + 1L
        variables$inf_num$queue_update(values = inf, index = target)
        # draw ab titre value
        inf[inf > length(parameters$mu_ab_infection)] <- length(parameters$mu_ab_infection)
        zdose <- log(10^rnorm(n = target$size(), mean = log10(parameters$mu_ab_infection[inf]),sd = std10_infection))
        zdose <- pmin(zdose, parameters$max_ab)
        variables$ab_titre$queue_update(values = zdose, index = target)
        # update last time of infection
        variables$inf_time$queue_update(values = timestep, index = target)
      }
    )
  }


}
