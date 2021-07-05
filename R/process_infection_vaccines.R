# --------------------------------------------------------------------------------
#   infection process for vaccination model (multiple doses, no types)
#   Sean L. Wu (slwood89@gmail.com)
#   July 2021
# --------------------------------------------------------------------------------

#' @title Infection process (nimue vaccine model)
#'
#' @description This samples infection events in the susceptible population. This incorporates
#' the slightly more complex force of infection calculation from the nimue model.
#'
#' @param parameters Model parameters
#' @param variables Model variable
#' @param events Model events
#' @param dt the time step
#' @export
infection_process_vaccine <- function(parameters, variables, events, dt) {

  stopifnot(all(c("states","vaccine_states","discrete_age") %in% names(variables)))

  doses_vec <- as.character(0:parameters$N_phase)
  doses_bset <- replicate(n = length(doses_vec),expr = {NULL})

  return(

    function(timestep) {

      infectious <- variables$states$get_index_of(c("IMild", "IAsymp", "ICase"))

      if (infectious$size() > 0) {

        day <- ceiling(timestep * dt)

        susceptible <- variables$states$get_index_of("S")

        # dose categories
        for (d in seq_along(doses_vec)) {
          doses_bset[[d]] <- variables$dose_num$get_index_if(doses_vec[d])
        }

        # age categories
        ages <- variables$discrete_age$get_values(infectious)

        # # infection by vaccine status
        # inf_vaxx <- variables$vaccine_states$get_values(index = infectious)
        #
        # # infection by age
        # ages <- variables$discrete_age$get_values(infectious)
        #
        # # compute cross tab for relative infectiousness, multiply by that matrix, and sum it out
        # inf_age_vax <- cross_tab_margins(a = ages,b = inf_vaxx,a_margin = 17,b_margin = 4)
        # inf_ages <- rowSums(inf_age_vax * parameters$rel_infectiousness_vaccinated)
        #
        # # calculate FoI for each age group
        # m <- get_contact_matrix(parameters)
        # lambda <- parameters$beta_set[day] * as.vector(m %*% (inf_ages * parameters$rel_infectiousness))
        #
        # # FoI for each susceptible person
        # sus_vaxx <- variables$vaccine_states$get_values(index = susceptible)
        # ages <- variables$discrete_age$get_values(susceptible)
        #
        # submat <- matrix(data = NA,nrow = susceptible$size(),ncol = 3)
        # submat[, 1] <- day
        # submat[, 2] <- ages
        # submat[, 3] <- sus_vaxx
        #
        # foi <- lambda[ages] * parameters$vaccine_efficacy_infection[submat]
        # susceptible$sample(rate = pexp(q = foi * dt))
        #
        # # newly infected persons queue the exposure event
        # if (susceptible$size() > 0) {
        #   events$exposure$schedule(susceptible, delay = 0)
        # }

      }

    } # end function

  ) # return
}
