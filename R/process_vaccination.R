# --------------------------------------------------
#   vaccination process for the general vaccine model
#   but with no types
#   Sean L. Wu (slwood89@gmail.com)
#   June 2021
# --------------------------------------------------

# on a day:
# 1. figure out what step we are on
# 2. allocate vaccines according to that step

vaccination_process <- function(parameters, variables, events, dt) {

  stopifnot(all(c("states","eligible","empty","dose_num","dose_time","phase","discrete_age") %in% names(variables)))

  return(

    function(timestep) {

      day <- ceiling(timestep * dt)

      # only distribute once a day
      if (day %% 1 == 0) {

        # get phase we are on
        phase <- variables$phase$value

        # how many doses we have today
        doses_left <- parameters$vaccine_set[day]

        # if vaccination done or no doses available, return early
        if (phase > parameters$N_phase || doses_left <= 0) {
          return(invisible(NULL))
        }

        # calculate what step we are on
        step <- get_vaccination_priority_stage(variables = variables, phase = phase, parameters = parameters)

        # advance to next phase: recalculate
        if (step == -1) {
          variables$phase$value <- variables$phase$value + 1L
          phase <- variables$phase$value
          # maybe don't need this? step should always be 1 at new phase
          step <- get_vaccination_priority_stage(variables = variables, phase = phase, parameters = parameters)
          stopifnot(step == 1)
          # need to do something when step = -1 and we are on the last phase; return early
          if (phase > parameters$N_phase) {
            return(invisible(NULL))
          }
        }

        stopifnot(phase <= parameters$N_phase)

        # row of the vaccine coverage matrix
        p_step <- vaccine_coverage_mat[step, ]

        targets <- target_pop(
          dose = phase, variables = variables, parameters = parameters, t = timestep, dt = dt,
          prioritisation = p_step,vaxx_priority = NULL
        )
        doses_given <- assign_doses(
          doses = doses_left,
          n_to_cover = targets$n_to_cover, eligible_age_bset = targets$eligible_age_bsets, eligible_age_counts = targets$eligible_age_counts,
          events = events,phase = phase,parameters = parameters
        )
        doses_left <- doses_left - doses_given

        # intermediate phases: give prioritized doses
        if (phase < parameters$N_phase) {

          targets_pri <- target_pop(
            dose = phase + 1, variables = variables, parameters = parameters, t = timestep, dt = dt,
            prioritisation = p_step,vaxx_priority = NULL
          )
          doses_today_pri <- assign_doses(
            doses = doses_left,
            n_to_cover = targets_pri$n_to_cover, eligible_age_bset = targets_pri$eligible_age_bsets, eligible_age_counts = targets_pri$eligible_age_counts,
            events = events,phase = phase + 1,parameters = parameters
          )
          doses_left <- doses_left - doses_given

        }

        # if remaining doses, give out for this dose according to the prioritization matrix
        if (doses_left > 0) {

          while( doses_left > 0 & step <= parameters$N_prioritisation_steps) {
            step <- step + 1
            p_step <- vaccine_coverage_mat[step, ]
            targets <- target_pop(
              dose = phase, variables = variables, parameters = parameters, t = timestep, dt = dt,
              prioritisation = p_step,vaxx_priority = NULL
            )
            doses_given <- assign_doses(
              doses = doses_left,
              n_to_cover = targets$n_to_cover, eligible_age_bset = targets$eligible_age_bsets, eligible_age_counts = targets$eligible_age_counts,
              events = events,phase = phase,parameters = parameters
            )
            doses_left <- doses_left - doses_given
          }

        } # end remaining doses check


      } # end vaccine distribution

    } # end function

  ) # end return

}
