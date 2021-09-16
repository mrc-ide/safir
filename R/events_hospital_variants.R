# --------------------------------------------------
#   hospitilisation scheduler listener with variants
#   Sean L. Wu (slwood89@gmail.com)
#   September 2021
# --------------------------------------------------


#' @title Create listener function to schedule events upon hospitilisation with variants of concern
#' @description
#' When the \code{hospitilisation} event fires, this listener should be called
#' to schedule future events and state changes for those persons.
#' @param parameters Model parameters
#' @param variables a list of all of the model variables
#' @param events a list of all of the model events
#' @export
create_hospital_scheduler_listener_voc <- function(
   parameters,
   variables,
   events
) {
  function(timestep, hospitalised) {

    disc_ages <- variables$discrete_age$get_values(hospitalised)
    prob_severe <- parameters$prob_severe[disc_ages]

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
            prob_death <- parameters$prob_severe_death_treatment[variables$discrete_age$get_values(mv_get)]
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
            prob_death <- parameters$prob_severe_death_no_treatment[variables$discrete_age$get_values(mv_not_get)]
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
            prob_death <- parameters$prob_non_severe_death_treatment[variables$discrete_age$get_values(ox_get)]
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
            prob_death <- parameters$prob_non_severe_death_no_treatment[variables$discrete_age$get_values(ox_not_get)]
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
