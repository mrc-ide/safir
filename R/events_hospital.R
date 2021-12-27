# --------------------------------------------------
#   hospitilisation scheduler listener (squire transmission model)
#   Sean L. Wu (slwood89@gmail.com)
#   May 2021
#   1. create_hospital_scheduler_listener
#   2. schedule_outcome
#   3. allocate_treatment
# --------------------------------------------------

# create_hospital_scheduler_listener_cpp <- function(
#   parameters,
#   variables,
#   events
# ) {
#   safir:::create_hospital_scheduler_listener_cpp_internal(
#     parameters = parameters,
#     states = variables$states$.variable,
#     discrete_age = variables$discrete_age$.variable,
#     imv_get_die = events$imv_get_die$.event,
#     imv_get_live = events$imv_get_live$.event,
#     imv_not_get_die = events$imv_not_get_die$.event,
#     imv_not_get_live = events$imv_not_get_live$.event,
#     iox_get_die = events$iox_get_die$.event,
#     iox_get_live = events$iox_get_live$.event,
#     iox_not_get_die = events$iox_not_get_die$.event,
#     iox_not_get_live = events$iox_not_get_live$.event
#   )
# }


#' @title Create listener function to schedule events upon hospitilisation
#' @description
#' When the \code{hospitilisation} event fires, this listener should be called
#' to schedule future events and state changes for those persons.
#' @param parameters Model parameters
#' @param variables a list of all of the model variables
#' @param events a list of all of the model events
#' @export
create_hospital_scheduler_listener <- function(
   parameters,
   variables,
   events
) {
  function(timestep, hospitalised) {

    day <- ceiling(timestep * parameters$dt)

    disc_ages <- variables$discrete_age$get_values(hospitalised)
    # prob_severe <- parameters$prob_severe[disc_ages]
    prob_severe <- get_probabilties(prob = parameters$prob_severe, ages = disc_ages, day = day)

    need_mv <- hospitalised$copy()
    need_mv$sample(prob_severe)

    need_ox <- hospitalised$set_difference(need_mv) # hospitalised should not be used after this point

    # individuals requiring mechanical ventilation (MV)
    if (need_mv$size() > 0) {

        # number of people who can get MV
        mv_get <- allocate_treatment(
            variables,
            need_treatment = need_mv,
            treated_state = c('IMVGetDie', 'IMVGetLive'),
            limit = parameters$ICU_beds
       )

        # schedule for those getting mv
        if (mv_get$size() > 0) {
            # prob_death <- parameters$prob_severe_death_treatment[variables$discrete_age$get_values(mv_get)]
            prob_death <- get_probabilties(prob = parameters$prob_severe_death_treatment, ages = variables$discrete_age$get_values(mv_get), day = day)
            schedule_outcome(
                target = mv_get,
                prob_successful = prob_death,
                success_event = events$imv_get_die,
                failure_event = events$imv_get_live
            )
        }

        # schedule for those not getting mv
        mv_not_get <- need_mv$set_difference(mv_get) # need_mv should not be used after this point
        if (mv_not_get$size() > 0) {
            # prob_death <- parameters$prob_severe_death_no_treatment[variables$discrete_age$get_values(mv_not_get)]
            prob_death <- get_probabilties(prob = parameters$prob_severe_death_no_treatment, ages = variables$discrete_age$get_values(mv_not_get), day = day)
            schedule_outcome(
                target = mv_not_get,
                prob_successful = prob_death,
                success_event = events$imv_not_get_die,
                failure_event = events$imv_not_get_live
            )
        }

    }

    # individuals requiring oxygen (Ox)
    if (need_ox$size() > 0) {

        # number of people who can get Ox
        ox_get <- allocate_treatment(
            variables,
            need_treatment = need_ox,
            treated_state = c('IOxGetDie', 'IOxGetLive', 'IRec'),
            limit = parameters$hosp_beds
        )

        # schedule for those getting ox
        if (ox_get$size() > 0) {
            # prob_death <- parameters$prob_non_severe_death_treatment[variables$discrete_age$get_values(ox_get)]
            prob_death <- get_probabilties(prob = parameters$prob_non_severe_death_treatment, ages = variables$discrete_age$get_values(ox_get), day = day)
            schedule_outcome(
                target = ox_get,
                prob_successful = prob_death,
                success_event = events$iox_get_die,
                failure_event = events$iox_get_live
            )
        }

        # schedule for those not getting ox
        ox_not_get <- need_ox$set_difference(ox_get) # need_ox should not be used after this point
        if (ox_not_get$size() > 0) {
            # prob_death <- parameters$prob_non_severe_death_no_treatment[variables$discrete_age$get_values(ox_not_get)]
            prob_death <- get_probabilties(prob = parameters$prob_non_severe_death_no_treatment, ages = variables$discrete_age$get_values(ox_not_get), day = day)
            schedule_outcome(
                target = ox_not_get,
                prob_successful = prob_death,
                success_event = events$iox_not_get_die,
                failure_event = events$iox_not_get_live
            )
        }

    }

    # end of function
  }
}


#' @title Schedule outcome
#' @description
#' schedule individuals into follow up events based based on bernoulli draws of
#' `prob_successful`
#' @param target the individuals to draw from
#' @param prob_successful the probability each target individual is successful
#' @param success_event will be scheduled on success
#' @param failure_event will be scheduled on failure
schedule_outcome <- function(
  target,
  prob_successful,
  success_event,
  failure_event
) {

    success <- target$copy()
    success$sample(prob_successful)

    failure <- target$copy()$set_difference(success)

    if (success$size() > 0) {
        success_event$schedule(target = success, delay = 0L)
    }

    if (failure$size() > 0) {
        failure_event$schedule(target = failure, delay = 0L)
    }

}


#' @title Allocate treatment
#' @description
#' sample a subset of individuals who will receive treatment. The subset is
#' allways smaller than the limit of treatments available minus those already
#' receiving treatment
#' @param variables Model variables
#' @param need_treatment a [individual::Bitset] of individuals who need treatment
#' @param treated_state a list states for individuals receiving treatment
#' @param limit the number of individuals who can receive treatment
#' @importFrom individual Bitset
allocate_treatment <- function(
  variables,
  need_treatment,
  treated_state,
  limit
) {
  # mv bed allocation
  occupied <- variables$states$get_index_of(treated_state)
  available <- limit - occupied$size()

  # who is getting an mv from available
  if (need_treatment$size() <= available) {
    return(need_treatment)
  }

  k <- max(0, available)
  if (k > 0) {
    get_treatment <- need_treatment$copy()
    get_treatment$choose(k = k)
    return(get_treatment)
  } else {
    return(Bitset$new(size = need_treatment$max_size))
  }
}


#' @noRd
get_probabilties <- function(prob, ages, day) {
  if (inherits(prob, "matrix")) {
    prob[cbind(ages, day)]
  } else {
    prob[ages]
  }
}
