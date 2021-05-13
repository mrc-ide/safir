#' @title Transitions from S
#'
#' @description S -> E if infected; can also become vaccinated.
#'
#' @param paramaters Model parameters
#' @param variables Model variable
#' @param events Model events
#' @noRd
S_process_zzz <- function(parameters, variables, events, vaccines = NULL, dt) {

    if (is.null(vaccines)) {

        function(timestep) {

            infectious <- variables$states$get_index_of(c("IMild", "IAsymp", "ICase"))

            if (infectious$size() > 0) {

                # Group infection by age
                ages <- variables$discrete_age$get_values(inf_states)
                inf_ages <- tabulate(ages, nbins = parameters$N_age)

                # calculate FoI for each age group
                m <- get_contact_matrix(parameters)
                lambda <- parameters$beta[as.integer(timestep * dt)] * rowSums(m %*% diag(inf_ages))

                # Transition from S to E
                susceptible <- variables$states$get_index_of("S")
                ages <- variables$discrete_age$get_values(susceptible)
                
                # FoI for each susceptible person
                lambda <- lambda[ages]

                # infected
                susceptible$sample(rate = pexp(q = lambda * dt))

                # newly infecteds queue the exposure event
                # no delay because of the discrete time model semantics; state must change at end of time step
                if (susceptible$size() > 0) {
                    events$exposure$schedule(susceptible, delay = 0)
                }
            }
        }

    } else {
        stop("not implemented yet!")
    }

}
