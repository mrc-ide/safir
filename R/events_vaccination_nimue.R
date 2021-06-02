# --------------------------------------------------
#   create event functions
#   May 2021
#   1. create_events
#   2. create_event_scheduler_listener
#   3. create_state_update_listener
#   4. attach_event_listeners
# --------------------------------------------------


#' @title Create events
#'
#' @param events a named list of individual::Event
#' @param parameters write me!
#' @param vaccines not currently used
#' @export
create_events_nimue <- function(events, parameters, vaccines = NULL) {

  # pop size
  N <- sum(parameters$population)

  events$v0_to_v1v2 = individual::TargetedEvent$new(N) # scheduled when vaccination occurs
  events$v1v2_to_v3v4 = individual::TargetedEvent$new(N) # scheduled when entering v1v2
  events$v3v4_to_v5 = individual::TargetedEvent$new(N) # scheduled when entering v3v4

  return(events)
}

# will need to attach additional listeners to recovery.
attach_event_listeners_nimue <- function(){

}
