# --------------------------------------------------------------------------------
#   vaccination process for nimue style vaccination model
#   Sean L. Wu (slwood89@gmail.com)
#   June 2021
# --------------------------------------------------------------------------------


#' @title Vaccination process (nimue vaccine model)
#'
#' @description This samples vaccination events (if there are vaccines available that day)
#' for individuals in S, E, R states.
#'
#' @param parameters Model parameters
#' @param variables Model variable
#' @param events Model events
#' @param dt the time step
#' @export
vaccination_process_nimue <- function(parameters, variables, events, dt) {

  stopifnot(all(c("states","eligible","vaccinated","empty","discrete_age") %in% names(variables)))
  stopifnot("v0_to_v1v2" %in% names(events))

  return(

    function(timestep) {

      mv <- parameters$vaccine_set[ceiling(timestep * dt)]

      if (mv > 0) {

        # calculate prioritisation step and which age groups are eligible right now
        pr <- vapply(X = 1:parameters$N_age, FUN = function(a){get_proportion_vaccinated_nimue(variables = variables,age = a)}, FUN.VALUE = numeric(1))

        vaccination_target_mat <- matrix(data = 0,nrow = parameters$N_prioritisation_steps,ncol = parameters$N_age)
        for (p in 1:parameters$N_prioritisation_steps) {
          vaccination_target_mat[p, ] <- as.integer(pr < parameters$vaccine_coverage_mat[p, ])
        }

        vaccine_target_vec <- rep(0, parameters$N_prioritisation_steps)
        for (p in 1:parameters$N_prioritisation_steps) {
          # an entire row summing to zero means that step has been completed
          vaccine_target_vec[p] <- as.integer(sum(vaccination_target_mat[p, ]) == 0)
        }
        current_index <- min(sum(vaccine_target_vec) + 1, parameters$N_prioritisation_steps)

        vaccination_target <- vaccination_target_mat[current_index, ]

        # if no vaccination targets remain don't run the code to distribute vaccines
        if (!all(vaccination_target == 0)) {

          # clear out eligible
          variables$eligible$and(variables$empty)

          SER <- variables$states$get_index_of(c("S","E","R"))
          target_ages <- which(vaccination_target > 0)
          SER$and(variables$discrete_age$get_index_of(target_ages))

          # set who is eligible: SER people in an age group in this priority step AND unvaccinated
          variables$eligible$or(SER)
          variables$eligible$and(variables$vaccinated$not())

          # calc rate of vaccination now
          vr_den <- variables$eligible$size()
          vr <- ifelse(mv <= 0, 0, min(mv/vr_den, 1))

          # sample who gets vaccinated
          variables$eligible$sample(rate = pexp(q = vr * dt))
          if (variables$eligible$size() > 0) {
            events$v0_to_v1v2$schedule(variables$eligible, delay = 0)
          }

        }

      } # end if

    } # end function

  ) # end return
}


#' @title C++ vaccination process (nimue vaccine model)
#'
#' @description This samples vaccination events (if there are vaccines available that day)
#' for individuals in S, E, R states.
#' Calls \code{\link{vaccination_process_nimue_cpp_internal}} to return an external pointer object.
#'
#' @param parameters Model parameters
#' @param variables Model variable
#' @param events Model events
#' @param dt the time step
#' @export
vaccination_process_nimue_cpp <- function(parameters, variables, events, dt) {

  stopifnot(all(c("states","eligible","vaccinated","empty","discrete_age") %in% names(variables)))
  stopifnot("v0_to_v1v2" %in% names(events))

  return(
    safir::vaccination_process_nimue_cpp_internal(
      parameters = parameters,
      states = variables$states$.variable,
      eligible = variables$eligible$.bitset,
      vaccinated = variables$vaccinated$.bitset,
      empty = variables$empty$.bitset,
      discrete_age = variables$discrete_age$.variable,
      v0_to_v1v2 = events$v0_to_v1v2$.event,
      dt = dt
    )
  )
}

