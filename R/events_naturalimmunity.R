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
#' @export
attach_event_listeners_natural_immunity <- function(variables, events, parameters, dt) {

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

  # boost antibody titre
  events$recovery$add_listener(
    function(timestep, target) {
      # update inf_num
      inf <- variables$inf_num$get_values(target) + 1L
      variables$inf_num$queue_update(values = inf, index = target)
      # draw ab titre value
      zdose <- log(10^rnorm(n = target$size(), mean = log10(parameters$mu_ab_infection[inf]),sd = parameters$std10))
      variables$ab_titre$queue_update(values = zdose, index = target)
      # update last time of infection
      variables$inf_time$queue_update(values = timestep, index = target)
    }
  )

}

# note that for Ab titre, the *only* parameter that is vaccine product specific
# is the boost upon each dose. So the decay rates are fixed (just protein decay).
# that's great.

# We really need to check that using variables$dose_time like this wont screw
# anything up. It will, because its used in distribution_vaccines.R to
# check if people are past the threshold. So we have to use something new.

# new stuff:
# 1. parameters$mu_ab_infection
# 2. variables$inf_num

# need to change:
# 1. vaccine_ab_titre_process: needs to update ab titre for people who are
# vaccinated OR infected. So do an AND with variables$inf_num > 0.

