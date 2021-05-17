# vaccination needs special handling, so all states from which individuals may be vaccinated get their own processes.


#' @title Transitions from S
#'
#' @description S -> E if infected; can also become vaccinated.
#'
#' @param paramaters Model parameters
#' @param variables Model variable
#' @param events Model events
#' @noRd
infection_process_zzz <- function(parameters, variables, events, dt, vaccines = NULL) {

  if (is.null(vaccines)) {

    return(

      # process without vaccination
      function(timestep) {

        infectious <- variables$states$get_index_of(c("IMild", "IAsymp", "ICase"))

        if (infectious$size() > 0) {

          # Group infection by age
          ages <- variables$discrete_age$get_values(infectious)
          inf_ages <- tabulate(ages, nbins = parameters$N_age)

          # calculate FoI for each age group
          m <- get_contact_matrix(parameters)
          lambda <- parameters$beta[1L + as.integer(timestep * dt)] * rowSums(m %*% diag(inf_ages))

          # Transition from S to E
          susceptible <- variables$states$get_index_of("S")
          ages <- variables$discrete_age$get_values(susceptible)

          # FoI for each susceptible person
          lambda <- lambda[ages]

          # infected
          susceptible$sample(rate = pexp(q = lambda * dt))

          # newly infecteds queue the exposure event
          if (susceptible$size() > 0) {
            events$exposure$schedule(susceptible, delay = 0)
          }

        }
      }

    )

  } else {
    stop("not implemented yet!")
  }

}

E_process_zzz <- function(parameters, variables, events, dt, shift = 1 , vaccines = NULL) {

  ICase_delay <- make_rerlang(mu = parameters$dur_E, dt = dt, shift = shift)
  IMild_delay <- make_rerlang(mu = parameters$dur_E, dt = dt, shift = shift)
  IAsymp_delay <- make_rerlang(mu = parameters$dur_E, dt = dt, shift = shift)

  # no vaccines
  if (is.null(vaccines)) {

    return(
      function(timestep) {

        exposed <- variables$states$get_index_of("E")

        ages <- variables$discrete_age$get_values(exposed)
        prob_hosp <- parameters$prob_hosp[ages]
        hosp <- bernoulli_multi_p(prob_hosp) # no exponential CDF because they go onto tracks immediately upon arrival in this state

        to_hosp <- which(hosp)
        no_hosp <- which(!hosp)

        # severe infections: E->ICase
        if (length(to_hosp) > 0) {
          events$severe_infection$schedule(
            target = individual::filter_bitset(exposed, to_hosp),
            delay = ICase_delay()
          )
        }

        if (length(no_hosp) > 0) {

          no_hosp_individuals <- individual::filter_bitset(exposed, no_hosp)

          # E->IAsymp
          prob_asymp <- parameters$prob_asymp[disc_ages[no_hosp]]
          asymp <- bernoulli_multi_p(prob_asymp)
          if (any(asymp)) {
            events$asymp_infection$schedule(
              target = individual::filter_bitset(exposed, which(asymp)),
              delay = IAsymp_delay()
            )
          }

          # E->IMild
          if (any(!asymp)) {
            events$asymp_infection$schedule(
              target = individual::filter_bitset(exposed, which(!asymp)),
              delay = IMild_delay()
            )
          }

        }

      }
    )

  # vaccines
  } else {
    stop("not implemented yet!")
  }

}

