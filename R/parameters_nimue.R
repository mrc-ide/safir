#' @title Append vaccine parameters from nimue model
#' @description This calls \code{\link[nimue]{parameters}} from \href{https://mrc-ide.github.io/nimue/index.html}{nimue}
#' and appends the following parameters to the list:
#'
#'  * gamma_vaccine_delay: mean duration of period from vaccination to vaccine protection.
#'  * gamma_V: mean duration of vaccine-derived immunity (days)
#'  * rel_infectiousness_vaccinated: Relative infectiousness per age category of vaccinated individuals relative to unvaccinated individuals. Default = rep(1, 17), which is no impact of vaccination on onwards transmissions
#'  * rel_infectiousness: Relative infectiousness per age category relative to maximum infectiousness category. Default = rep(1, 17)
#'  * vaccine_efficacy_infection: Efficacy of vaccine against infection. This parameter must either be length 1 numeric (a single efficacy for all age groups) or length 17 numeric vector (an efficacy for each age group). An efficacy of 1 will reduce FOI by 100 percent, an efficacy of 0.2 will reduce FOI by 20 percent etc.
#'  * prob_hosp: Efficacy of vaccine against severe (requiring hospitilisation) disease (by age). This parameter must either be length 1 numeric (a single efficacy for all age groups) or length 17 numeric vector (an efficacy for each age group). An efficacy of 1 will reduce the probability of hospitalisation by 100 percent, an efficacy of 0.2 will reduce the probability of hospitalisation by 20 percent etc.
#'  * vaccine_coverage_mat: Vaccine coverage targets by age (columns) and priority (row)
#'  * N_prioritisation_steps: number of vaccine prioritization levels
#'  * current_prioritisation_step: current prioritization step we are on
#'  * vaccine_set: vaccines available each day of simulation
#'
#' @param parameters list from [get_parameters]
#' @param ... Other parameters for \code{nimue:::}\code{\link[nimue]{parameters}}
#' @export
append_vaccine_nimue <- function(parameters, ...) {

  nimue_pars <- call_nimue_pars(country = parameters$country, ...)

  parameters$gamma_vaccine_delay <- 1 / (nimue_pars$gamma_vaccine[2] / 2)
  parameters$gamma_V <- 1 / (nimue_pars$gamma_vaccine[4] / 2)
  parameters$rel_infectiousness_vaccinated <- nimue_pars$rel_infectiousness_vaccinated[,-c(3,5)]
  parameters$rel_infectiousness <- nimue_pars$rel_infectiousness
  parameters$vaccine_coverage_mat <- nimue_pars$vaccine_coverage_mat
  parameters$N_prioritisation_steps <- nimue_pars$N_prioritisation_steps
  parameters$current_prioritisation_step <- 1L

  parameters$vaccine_set <- interp_input_par(
    x = c(nimue_pars$tt_vaccine, parameters$time_period),
    y = c(nimue_pars$max_vaccine, tail(nimue_pars$max_vaccine, 1))
  )

  parameters$prob_hosp <- interp_nimue_array(
    x = c(nimue_pars$tt_vaccine_efficacy_disease, parameters$time_period),
    y = nimue_pars$prob_hosp
  )

  parameters$vaccine_efficacy_infection <- interp_nimue_array(
    x = c(nimue_pars$tt_vaccine_efficacy_infection, parameters$time_period),
    y = nimue_pars$vaccine_efficacy_infection
  )

  parameters$prob_severe_death_no_treatment <- nimue_pars$prob_severe_death_no_treatment
  parameters$prob_severe_death_treatment <- nimue_pars$prob_severe_death_treatment
  parameters$prob_non_severe_death_no_treatment <- nimue_pars$prob_non_severe_death_no_treatment
  parameters$prob_non_severe_death_treatment <- nimue_pars$prob_non_severe_death_treatment

  return(parameters)
}


#' @noRd
call_nimue_pars <- function(
  country,
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
  nimue:::parameters(
    country = country,
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
