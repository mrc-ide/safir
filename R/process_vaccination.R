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

        # calculate what step we are on
        step <- get_vaccination_priority_stage(variables = variables, phase = phase, parameters = parameters)

        # advance to next phase: recalculate
        if (step == -1) {
          variables$phase$value <- variables$phase$value + 1L
          phase <- variables$phase$value
          # maybe don't need this? step should always be 1 at new phase
          step <- get_vaccination_priority_stage(variables = variables, phase = phase, parameters = parameters)
          # need to do something when step = -1 and we are on the last phase
        }

        stopifnot(phase <= parameters$N_phase)

        # how many doses we have today
        doses_today <- parameters$vaccine_set[day]

        p_step <- vaccine_coverage_mat[step, ]

        # intermediate phases: give prioritized doses
        if (phase < parameters$N_phase) {

          targets <- target_pop(
            dose = phase, variables = variables, parameters = parameters, t = timestep, dt = dt, prioritisation = p_step,vaxx_priority = NULL
          )
          safir::assign_doses(
            t = t,dt = 1,doses = doses,
            n_to_cover = targeted$n_to_cover,eligible_age_bset = targeted$eligible_age_bsets,eligible_age_counts = targeted$eligible_age_counts,
            variables = variables,events = events,phase = phase,parameters = parameters
          )

        # final phase: do not give prioritized doses for phase + 1
        } else {

        }



      } # end vaccine distribution





    } # end function

  ) # end return

}
