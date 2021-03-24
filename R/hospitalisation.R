#' @title Simulate hospitilisation flow
#' @description
#' Allocated individuals to hospital resources upon arrival
#' @param hospitalised vector of individuals who need hospitalisation
#' @param discrete_age the discrete age variable
#' @param human humans
#' @param states a list of all of the model states
#' @param events a list of all of the model events
#' @noRd
hospitilisation_flow_process <- function(
   discrete_age,
   human,
   states,
   events
) {
  function(api, hospitalised) {
    parameters <- api$get_parameters()
    disc_ages <- api$get_variable(human, discrete_age)
    prob_severe <- prob_outcome(hospitalised, disc_ages,
                                parameters$prob_severe)

    # 1. Who needs a MV
    mv_success <- bernoulli_multi_p(prob_severe)
    need_mv <- hospitalised[mv_success]
    if(length(need_mv) > 0) {
       mv_get <- allocate_treatment(
         api = api,
         human = human,
         need_treatment = need_mv,
         treated_state = states[c('IMVGetDie', 'IMVGetLive')],
         limit = parameters$ICU_beds
       )

       # schedule for those getting mv
       if (length(mv_get) > 0) {
         schedule_outcome(
           api = api,
           target = mv_get,
           prob_successful = prob_outcome(
             mv_get,
             disc_ages,
             parameters$prob_severe_death_treatment
           ),
           success_event = events$imv_get_die,
           failure_event = events$imv_get_live
         )
       }

       # schedule for those not getting mv
       mv_not_get <- setdiff(need_mv, mv_get)
       if (length(mv_not_get) > 0) {
         schedule_outcome(
           api = api,
           target = mv_not_get,
           prob_successful = prob_outcome(
             mv_not_get,
             disc_ages,
             parameters$prob_severe_death_no_treatment
           ),
           success_event = events$imv_not_get_die,
           failure_event = events$imv_not_get_live
         )
       }
    }

    # 2. Who needs Ox
    need_ox <- hospitalised[!mv_success]
    if (length(need_ox) > 0) {

      ox_get <- allocate_treatment(
        api = api,
        human = human,
        need_treatment = need_ox,
        treated_state = states[c('IOxGetDie', 'IOxGetLive', 'IRec')],
        limit = parameters$hosp_beds
      )

      # schedule for those getting ox
      if (length(ox_get) > 0) {
        schedule_outcome(
          api = api,
          target = ox_get,
          prob_successful = prob_outcome(
            ox_get,
            disc_ages,
            parameters$prob_non_severe_death_treatment
          ),
          success_event = events$iox_get_die,
          failure_event = events$iox_get_live
        )
      }

      # schedule for those not getting ox
      ox_not_get <- setdiff(need_ox, ox_get)
      if (length(ox_not_get) > 0) {
        schedule_outcome(
          api = api,
          target = ox_not_get,
          prob_successful = prob_outcome(
            ox_not_get,
            disc_ages,
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
#' @param api simulation api
#' @param human humans
#' @param need_treatment a vector of individuals who need treatment
#' @param treated_state a list states for individuals receiving treatment
#' @param limit the number of individuals who can receive treatment
#' @noRd
allocate_treatment <- function(
  api,
  human,
  need_treatment,
  treated_state,
  limit
) {
  # mv bed allocation
  occupied <- api$get_state(human, treated_state)
  available <- limit - length(occupied)

  # who is getting an mv from available
  if (length(need_treatment) <= available) {
    return(need_treatment)
  }

  need_treatment[sample.int(length(need_treatment), available)]
}

#' @title Schedule outcome
#' @description
#' schedule individuals into follow up events based based on bernoulli draws of
#' `prob_successful`
#' @param api simulation api
#' @param target the individuals to draw from
#' @param prob_successful the probability each target individual is successful
#' @param success_event will be scheduled on success
#' @param failure_event will be scheduled on failure
#' @noRd
schedule_outcome <- function(
  api,
  target,
  prob_successful,
  success_event,
  failure_event
) {
  success <- bernoulli_multi_p(prob_successful)

  if(sum(success) > 0) {
    api$schedule(
      event = success_event,
      target = target[success],
      delay = 0
    )
  }

  if(sum(!success) > 0) {
    api$schedule(
      event = failure_event,
      target = target[!success],
      delay = 0
    )
  }
}

#' @title Probability of outcome
#' @description
#' get the probabilities of an outcome for a target population
#' @param target the individuals of interest
#' @param age the discrete age group of all individuals
#' @param probs the probabilities per age group
#' @noRd
prob_outcome <- function(
  target,
  age,
  probs
) {
  probs[as.integer(age[target])]
}
