# --------------------------------------------------
#   vaccination process for the general vaccine model
#   but with no types
#   Sean L. Wu (slwood89@gmail.com)
#   June 2021
# --------------------------------------------------

# on a day:
# 1. figure out what step we are on
# 2. allocate vaccines according to that step

#' @title Vaccination distribution process (multiple doses, no types)
#' @description This function distributes vaccine doses each day for the vaccine
#' model with multiple doses but no types. It distributes vaccines according to the
#' following strategy.
#'
#'     * start at \code{phase = 1}, meaning we begin by distributing the first dose
#'     * use \code{\link{get_vaccination_priority_stage}} to check what stage of the
#'     prioritization matrix \code{vaccine_coverage_mat} from \code{\link[nimue]{strategy_matrix}}
#'     we should be using to decide who to vaccinate.
#'     * we can advance one row of the prioritization matrix when coverage for this dose (phase) >=
#'     prioritization matrix and coverage for next dose (phase) >= prioritization matrix groups prioritized
#'     for the next dose by the matrix \code{next_dose_priority}.
#'
#' @param parameters a named list
#' @param variables a named list
#' @param events a named list
#' @param dt size of time step
#' @export
vaccination_process <- function(parameters, variables, events, dt) {

  stopifnot(all(c("dose_num","dose_time","phase","discrete_age") %in% names(variables)))
  stopifnot(all(c("N_phase", "dose_period", "vaccine_coverage_mat", "N_prioritisation_steps", "next_dose_priority", "vaccine_set") %in% names(parameters)))

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
        p_step <- parameters$vaccine_coverage_mat[step, ]

        # assign these doses to the current dose phase and priority group
        targets <- target_pop(
          dose = phase, variables = variables, parameters = parameters, t = timestep, dt = dt,
          prioritisation = p_step,vaxx_priority = NULL
        )
        doses_left <- assign_doses(
          doses = doses_left,
          n_to_cover = targets$n_to_cover, eligible_age_bset = targets$eligible_age_bsets, eligible_age_counts = targets$eligible_age_counts,
          events = events,dose = phase,parameters = parameters
        )

        # intermediate phases: give prioritized doses to next group for next phase
        if (phase < parameters$N_phase & doses_left > 0) {

          targets_pri <- target_pop(
            dose = phase + 1, variables = variables, parameters = parameters, t = timestep, dt = dt,
            prioritisation = p_step,vaxx_priority = parameters$next_dose_priority[phase, ]
          )
          doses_left <- assign_doses(
            doses = doses_left,
            n_to_cover = targets_pri$n_to_cover, eligible_age_bset = targets_pri$eligible_age_bsets, eligible_age_counts = targets_pri$eligible_age_counts,
            events = events,dose = phase + 1,parameters = parameters
          )

        }

        # if remaining doses, give out for this dose phase according to the prioritization matrix
        if (doses_left > 0) {

          while( doses_left > 0 & step <= parameters$N_prioritisation_steps) {
            step <- step + 1
            p_step <- parameters$vaccine_coverage_mat[step, ]
            targets <- target_pop(
              dose = phase, variables = variables, parameters = parameters, t = timestep, dt = dt,
              prioritisation = p_step,vaxx_priority = NULL
            )
            doses_left <- assign_doses(
              doses = doses_left,
              n_to_cover = targets$n_to_cover, eligible_age_bset = targets$eligible_age_bsets, eligible_age_counts = targets$eligible_age_counts,
              events = events,dose = phase,parameters = parameters
            )
          }

        } # end remaining doses check

      } # end vaccine distribution

    } # end function

  ) # end return

}
