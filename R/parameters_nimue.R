#  --------------------------------------------------------------------------------
#   get parameters for nimue vaccination model
#   Sean L. Wu (slwood89@gmail.com)
#   June 2021
#  --------------------------------------------------------------------------------

#' @title Get parameters for safir (nimue vaccine model)
#' @description This calls \code{\link[nimue]{parameters}} from \href{https://mrc-ide.github.io/nimue/index.html}{nimue}.
#'
#' @param population integer vector of population by age group
#' @param contact_mat a contact matrix from \code{\link[squire]{get_mixing_matrix}}
#' @param time_period maximum time of simulation
#' @param max_age max age for humans
#' @param lambda_external a vector of length equal to \code{time_period} giving an additional additive term contributing to the force of infection.
#' Models infectious contacts the population has with external (unmodeled) populations.
#' @param dt time step size
#' @param ... Other parameters for [nimue::parameters]
#' @export
get_parameters_nimue <- function(
  population,
  contact_mat,
  time_period,
  max_age = 100,
  lambda_external = NULL,
  dt = 1,
  ...
) {

  stopifnot(length(population)==17)
  stopifnot(is.finite(population))

  pars <- call_nimue_pars(
    population = population,
    contact_matrix_set = contact_mat,
    time_period = time_period,
    ...
  )

  pars$dt <- dt
  pars$max_age <- max_age
  pars$time_period <- time_period

  # time varying

  # nimue vaccine supply
  pars$vaccine_set <- interp_input_par(
    x = c(pars$tt_vaccine, time_period),
    y = c(pars$max_vaccine, tail(pars$max_vaccine, 1))
  )

  # nimue probability of hospitalization (from vaccine_efficacy_disease)
  pars$prob_hosp <- interp_nimue_array(
    x = c(pars$tt_vaccine_efficacy_disease, time_period),
    y = pars$prob_hosp
  )

  # nimue probability of infection via vaxx
  pars$vaccine_efficacy_infection <- interp_nimue_array(
    x = c(pars$tt_vaccine_efficacy_infection, time_period),
    y = pars$vaccine_efficacy_infection
  )

  # beta
  pars$beta_set <- interp_input_par(
    c(pars$tt_beta, time_period),
    c(pars$beta_set, tail(pars$beta_set, 1))
  )

  # default asymptomatic pars
  asymp_pars <- get_asymptomatic()
  pars$dur_IAsymp <- asymp_pars$dur_IAsymp
  pars$prob_asymp <- asymp_pars$prob_asymp
  pars$IAsymp_0 <- asymp_pars$IAsymp_0

  # turn all rates back into durations
  pars$dur_E <- 2 / pars$gamma_E
  pars$dur_IMild <- 1 / pars$gamma_IMild
  pars$dur_ICase <- 2 / pars$gamma_ICase
  pars$dur_get_ox_survive <- 2 / pars$gamma_get_ox_survive
  pars$dur_get_ox_die <- 2 / pars$gamma_get_ox_die
  pars$dur_not_get_ox_survive <- 2 / pars$gamma_not_get_ox_survive
  pars$dur_not_get_ox_die <- 2 / pars$gamma_not_get_ox_die
  pars$dur_get_mv_survive <- 2 / pars$gamma_get_mv_survive
  pars$dur_get_mv_die <- 2 / pars$gamma_get_mv_die
  pars$dur_not_get_mv_survive <- 2 / pars$gamma_not_get_mv_survive
  pars$dur_not_get_mv_die <- 2 / pars$gamma_not_get_mv_die
  pars$dur_rec <- 2 / pars$gamma_rec

  pars$dur_R <- 2 / pars$gamma_R
  pars$dur_V <- 2 / pars$gamma_vaccine[4]
  pars$dur_vaccine_delay <- 2 / pars$gamma_vaccine[2]

  # drop extra cols to line up with 1:4 vaccine states
  pars$rel_infectiousness_vaccinated <- pars$rel_infectiousness_vaccinated[,-c(3,5)]

  # external FoI term
  if (!is.null(lambda_external)) {
    stopifnot(length(lambda_external) == time_period)
    stopifnot(all(is.finite(lambda_external)))
    stopifnot(all(lambda_external >= 0))
    pars$lambda_external <- as.numeric(lambda_external)
  } else {
    pars$lambda_external <- rep(0, time_period)
  }

  return(pars)
}


#' @noRd
call_nimue_pars <- function(
  population,
  contact_matrix_set,
  # transmission
  R0 = 3,
  tt_R0 = 0,
  beta_set = NULL,
  # initial state, duration, reps
  time_period,
  # probabilities
  prob_hosp = nimue:::probs$prob_hosp,
  prob_severe = nimue:::probs$prob_severe,
  prob_non_severe_death_treatment = nimue:::probs$prob_non_severe_death_treatment,
  prob_non_severe_death_no_treatment = nimue:::probs$prob_non_severe_death_no_treatment,
  prob_severe_death_treatment = nimue:::probs$prob_severe_death_treatment,
  prob_severe_death_no_treatment = nimue:::probs$prob_severe_death_no_treatment,
  p_dist = nimue:::probs$p_dist,
  # onward infectiousness
  rel_infectiousness = nimue:::probs$rel_infectiousness,
  rel_infectiousness_vaccinated = nimue:::probs$rel_infectiousness_vaccinated,
  # durations
  dur_E  = nimue:::durs$dur_E,
  dur_IMild = nimue:::durs$dur_IMild,
  dur_ICase = nimue:::durs$dur_ICase,
  # hospital durations
  dur_get_ox_survive = nimue:::durs$dur_get_ox_survive,
  dur_get_ox_die = nimue:::durs$dur_get_ox_die,
  dur_not_get_ox_survive = nimue:::durs$dur_not_get_ox_survive,
  dur_not_get_ox_die = nimue:::durs$dur_not_get_ox_die,
  dur_get_mv_survive = nimue:::durs$dur_get_mv_survive,
  dur_get_mv_die = nimue:::durs$dur_get_mv_die,
  dur_not_get_mv_survive = nimue:::durs$dur_not_get_mv_survive,
  dur_not_get_mv_die = nimue:::durs$dur_not_get_mv_die,
  dur_rec = nimue:::durs$dur_rec,
  # vaccine
  dur_R = nimue:::vaccine_pars$dur_R,
  tt_dur_R = nimue:::vaccine_pars$tt_dur_R,
  dur_V = nimue:::vaccine_pars$dur_V,
  vaccine_efficacy_infection = nimue:::vaccine_pars$vaccine_efficacy_infection,
  tt_vaccine_efficacy_infection = nimue:::vaccine_pars$tt_vaccine_efficacy_infection,
  vaccine_efficacy_disease = nimue:::vaccine_pars$vaccine_efficacy_disease,
  tt_vaccine_efficacy_disease = nimue:::vaccine_pars$tt_vaccine_efficacy_disease,
  max_vaccine = nimue:::vaccine_pars$max_vaccine,
  tt_vaccine = nimue:::vaccine_pars$tt_vaccine,
  dur_vaccine_delay = nimue:::vaccine_pars$dur_vaccine_delay,
  vaccine_coverage_mat = nimue:::vaccine_pars$vaccine_coverage_mat,
  # health system capacity
  hosp_bed_capacity = NULL,
  ICU_bed_capacity = NULL,
  tt_hosp_beds = 0,
  tt_ICU_beds = 0,
  # misc
  seeding_cases = 20,
  seeding_age_order = NULL
) {
  nimue_parameters(
    country = NULL,
    population = population,
    contact_matrix_set = contact_matrix_set,
    # transmission
    R0 = R0,
    tt_R0 = tt_R0,
    beta_set = beta_set,
    # probabilities
    prob_hosp = prob_hosp,
    prob_severe = prob_severe,
    prob_non_severe_death_treatment = prob_non_severe_death_treatment,
    prob_non_severe_death_no_treatment = prob_non_severe_death_no_treatment,
    prob_severe_death_treatment = prob_severe_death_treatment,
    prob_severe_death_no_treatment = prob_severe_death_no_treatment,
    p_dist = p_dist,
    # onward infectiousness
    rel_infectiousness = rel_infectiousness,
    rel_infectiousness_vaccinated = rel_infectiousness_vaccinated,
    # durations
    dur_E  = dur_E,
    dur_IMild = dur_IMild,
    dur_ICase = dur_ICase,
    # hospital durations
    dur_get_ox_survive = dur_get_ox_survive,
    dur_get_ox_die = dur_get_ox_die,
    dur_not_get_ox_survive = dur_not_get_ox_survive,
    dur_not_get_ox_die = dur_not_get_ox_die,
    dur_get_mv_survive = dur_get_mv_survive,
    dur_get_mv_die = dur_get_mv_die,
    dur_not_get_mv_survive = dur_not_get_mv_survive,
    dur_not_get_mv_die = dur_not_get_mv_die,
    dur_rec = dur_rec,
    # vaccine
    dur_R = dur_R,
    tt_dur_R = tt_dur_R,
    dur_V = dur_V,
    vaccine_efficacy_infection = vaccine_efficacy_infection,
    tt_vaccine_efficacy_infection = tt_vaccine_efficacy_infection,
    vaccine_efficacy_disease = vaccine_efficacy_disease,
    tt_vaccine_efficacy_disease = tt_vaccine_efficacy_disease,
    max_vaccine = max_vaccine,
    tt_vaccine = tt_vaccine,
    dur_vaccine_delay = dur_vaccine_delay,
    vaccine_coverage_mat = vaccine_coverage_mat,
    # health system capacity
    hosp_bed_capacity = hosp_bed_capacity,
    ICU_bed_capacity = ICU_bed_capacity,
    tt_hosp_beds = tt_hosp_beds,
    tt_ICU_beds = tt_ICU_beds,
    # misc
    seeding_cases = seeding_cases,
    seeding_age_order = seeding_age_order
  )
}


#' @noRd
interp_nimue_array <- function(x, y, by = 1, end = max(x) + 1) {

  stopifnot(inherits(y, "array"))
  stopifnot(length(dim(y)) == 3)
  stopifnot(dim(y)[2] == 17)
  stopifnot(dim(y)[3] == 6)
  stopifnot(dim(y)[1] >= 1)

  # create the time series
  if (length(x) > 1) {
    x_all <- c(x[-length(x)], end)
  } else {
    x_all <- c(x, end)
  }

  xout <- seq(1, end, 1)
  yout <- array(data = NA,dim = c(length(xout), 17, 4))

  xdiffs <- diff(x_all)

  vacc_stages <- c(1, 2, 4, 6)
  for (i in seq_along(vacc_stages)) {
    for (j in seq_along(xdiffs)) {
      yout[(x_all+1)[j]:x_all[j+1], , i] <- do.call(rbind,replicate(n = xdiffs[j],expr = y[j,,vacc_stages[i]],simplify = FALSE))
    }
  }

  return(yout)
}
