#' @title Simulate hospitilisation flow
#' @description
#' Allocated individuals to hospital resources upon arrival
#' @param parameters Model parameters
#' @param states a list of all of the model states
#' @param events a list of all of the model events
#' @noRd
hospitilisation_flow_process <- function(
   parameters,
   variables,
   events
) {
  function(timestep, hospitalised) {
    disc_ages <- variables$discrete_age$get_values(hospitalised)
    prob_severe <- prob_outcome(disc_ages,
                                parameters$prob_severe)

    # 1. Who needs a MV
    mv_success <- bernoulli_multi_p(prob_severe)
    need_mv <- individual::filter_bitset(hospitalised, which(mv_success))
    if(need_mv$size() > 0) {
      browser()
       mv_get <- allocate_treatment(
         variables,
         need_treatment = need_mv,
         treated_state = c('IMVGetDie', 'IMVGetLive'),
         limit = parameters$ICU_beds
       )

       # schedule for those getting mv
       if (mv_get$size() > 0) {
         schedule_outcome(
           target = mv_get,
           prob_successful = prob_outcome(
             variables$discrete_age$get_values(mv_get),
             parameters$prob_severe_death_treatment
           ),
           success_event = events$imv_get_die,
           failure_event = events$imv_get_live
         )
       }

       # schedule for those not getting mv
       mv_not_get <- need_mv$set_difference(mv_get)
       if (length(mv_not_get) > 0) {
         schedule_outcome(
           target = mv_not_get,
           prob_successful = prob_outcome(
             variables$discrete_age$get_values(mv_not_get),
             parameters$prob_severe_death_no_treatment
           ),
           success_event = events$imv_not_get_die,
           failure_event = events$imv_not_get_live
         )
       }
    }

    # 2. Who needs Ox
    need_ox <- individual::filter_bitset(hospitalised, which(!mv_success))
    if (need_ox$size() > 0) {

      ox_get <- allocate_treatment(
        variables,
        need_treatment = need_ox,
        treated_state = c('IOxGetDie', 'IOxGetLive', 'IRec'),
        limit = parameters$hosp_beds
      )

      # schedule for those getting ox
      if (ox_get$size() > 0) {
        schedule_outcome(
          target = ox_get,
          prob_successful = prob_outcome(
            variables$discrete_age$get_values(ox_get),
            parameters$prob_non_severe_death_treatment
          ),
          success_event = events$iox_get_die,
          failure_event = events$iox_get_live
        )
      }
      # schedule for those not getting ox
      ox_not_get <- need_ox$set_difference(ox_get)
      if (ox_not_get$size() > 0) {
        schedule_outcome(
          target = ox_not_get,
          prob_successful = prob_outcome(
            variables$discrete_age$get_values(ox_not_get),
            parameters$prob_non_severe_death_no_treatment
          ),
          success_event = events$iox_not_get_die,
          failure_event = events$iox_not_get_live
        )
      }
    }
  }
}

#' @title Allocate treatment
#' @description
#' sample a subset of individuals who will receive treatment. The subset is
#' allways smaller than the limit of treatments available minus those already
#' receiving treatment
#' @param variables Model variables
#' @param need_treatment a vector of individuals who need treatment
#' @param treated_state a list states for individuals receiving treatment
#' @param limit the number of individuals who can receive treatment
#' @noRd
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

  get_treatment <- individual::filter_bitset(need_treatment, sample.int(need_treatment$size(), max(0, available)))
  return(get_treatment)
}

#' @title Schedule outcome
#' @description
#' schedule individuals into follow up events based based on bernoulli draws of
#' `prob_successful`
#' @param target the individuals to draw from
#' @param prob_successful the probability each target individual is successful
#' @param success_event will be scheduled on success
#' @param failure_event will be scheduled on failure
#' @noRd
schedule_outcome <- function(
  target,
  prob_successful,
  success_event,
  failure_event
) {
  success <- bernoulli_multi_p(prob_successful)

  if(sum(success) > 0) {
    to_success <- individual::filter_bitset(target, which(success))
    success_event$schedule(to_success, delay = 0)
  }

  if(sum(!success) > 0) {
    to_failure <- individual::filter_bitset(target, which(!success))
    failure_event$schedule(to_failure, delay = 0)
  }
}

#' @title Probability of outcome
#' @description
#' get the probabilities of an outcome for a target population
#' @param age the discrete age group of all individuals
#' @param probs the probabilities per age group
#' @noRd
prob_outcome <- function(
  age,
  probs
) {
  probs[age]
}
