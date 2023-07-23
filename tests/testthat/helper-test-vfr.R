draw_nt_vfr <- function(parameters, n, tmax, vfr, vfr_time_1, vfr_time_2) {

  # get pars out
  hl_s <- parameters$hl_s # Half life of antibody decay - short
  hl_l <- parameters$hl_l # Half life of antibody decay - long
  period_s <- parameters$period_s
  t_period_l <- parameters$t_period_l # Time point at which to have switched to longest half-life

  time_to_decay <- t_period_l - period_s # time in days to reach longest half-life
  dr_s <- -log(2)/hl_s # Corresponding decay rate in days for half life above
  dr_l <- -log(2)/hl_l

  ab_50 <- parameters$ab_50 # titre relative to convalescent required to provide 50% protection from infection, on linear scale
  ab_50_severe <- parameters$ab_50_severe
  k <- parameters$k # shape parameter of efficacy curve

  mu_ab_d1 <- parameters$mu_ab[1] # mean titre dose 1
  std10 <- parameters$std10 # Pooled standard deviation of antibody level on log10 data

  t <- 0:tmax

  # vector of decay rates over time: first we have a period of fast decay, then gradually shift to period of long decay
  dr_vec <- c(rep(dr_s, period_s),
              seq(dr_s, dr_l, length.out = time_to_decay),
              rep(dr_l, (length(t) - t_period_l)))

  z1 <- rnorm(n, log10(mu_ab_d1), std10)

  # initiate titre vector
  nt <- matrix(data = NaN,nrow = length(t),ncol = n)
  nt[1, ] <- log(10^z1)

  # decay antibodies over time on natural log scale
  for (i in (2:length(t))){
    nt[i, ] <- nt[i-1, ] + dr_vec[i-1]
  }

  nt_log <- nt

  # VFR
  vfr_vector <- c(rep(1, times = vfr_time_1), seq(from = 1, to = vfr, length.out = (vfr_time_2 - vfr_time_1 + 1)), rep(vfr, times = tmax - vfr_time_2))

  nt <- exp(nt) # return to linear scale
  nt <- nt / vfr_vector

  # relate titre to efficacy over time - using log-10 parameters
  ef_infection <- 1 / (1 + exp(-k * (log10(nt) - log10(ab_50))))
  ef_severe_uncond <- 1 / (1 + exp(-k * (log10(nt) - log10(ab_50_severe))))
  ef_severe <-  1 - ((1 - ef_severe_uncond)/(1 - ef_infection))

  return(list(nt = nt_log, z1 = z1, ef_infection = ef_infection, ef_severe = ef_severe))
}



simulate_vfr <- function(iso3c, vfr, tmax, dt, R0, ab_titre, pop, mu_ab_infection = NULL, ret_ab = FALSE) {

  contact_mat <- squire::get_mixing_matrix(iso3c = iso3c)

  # vaccine dosing
  vaccine_doses <- 2
  dose_period <- c(NaN, 28)
  vaccine_set <- rep(0, tmax)

  # vaccine strategy
  vaccine_coverage_mat <- nimue::strategy_matrix(strategy = "Elderly",max_coverage = 0.2)
  next_dose_priority <- matrix(data = 0, nrow = vaccine_doses - 1,ncol = ncol(vaccine_coverage_mat))
  next_dose_priority[1, 15:17] <- 1 # prioritize 3 oldest age groups for next dose

  # base parameters
  parameters <- safir::get_parameters(
    population = pop$n,
    contact_matrix_set = contact_mat,
    iso3c = iso3c,
    R0 = R0,
    time_period = tmax,
    dt = dt
  )

  # vaccine parameters
  ab_parameters <- get_vaccine_ab_titre_parameters(vaccine = "Pfizer", max_dose = vaccine_doses,correlated = FALSE)

  # combine parameters and verify
  parameters <- make_vaccine_parameters(
    safir_parameters = parameters,
    vaccine_ab_parameters = ab_parameters,
    vaccine_set = vaccine_set,
    dose_period = dose_period,
    strategy_matrix = vaccine_coverage_mat,
    next_dose_priority_matrix = next_dose_priority
  )

  parameters <- make_immune_parameters(parameters = parameters, vfr = vfr, mu_ab_infection = mu_ab_infection)

  # create variables
  timesteps <- parameters$time_period/dt
  variables <- create_variables(pop = pop, parameters = parameters)
  variables <- create_vaccine_variables(variables = variables,parameters = parameters)
  variables <- create_natural_immunity_variables(variables = variables, parameters = parameters)
  variables <- create_independent_nat_variables(variables = variables, parameters = parameters)

  # create events
  events <- create_events(parameters = parameters)
  events <- create_events_vaccination(events = events,parameters = parameters)
  attach_event_listeners(variables = variables,events = events,parameters = parameters, dt = dt)
  attach_event_listeners_vaccination(variables = variables,events = events,parameters = parameters,dt = dt)
  attach_event_listeners_natural_immunity(variables = variables, events = events, parameters = parameters, dt = dt, additive = T)

  # make renderers
  renderer <- individual::Render$new(parameters$time_period)

  # processes
  processes <- list(
    natural_immunity_ab_titre_process(parameters = parameters,variables = variables,dt = dt),
    vaccination_process(parameters = parameters,variables = variables,events = events,dt = dt),
    infection_process_vaccine_cpp(parameters = parameters,variables = variables,events = events,dt = dt),
    categorical_count_renderer_process_daily(renderer = renderer,variable = variables$states,categories = variables$states$get_categories(),dt = dt)
  )

  setup_events(parameters = parameters,events = events,variables = variables,dt = dt)

  # give everyone ab titre
  variables$ab_titre$queue_update(values = ab_titre, index = 1:sum(pop$n))
  variables$inf_time$queue_update(values = 0, index = 1:sum(pop$n))
  variables$inf_num$queue_update(values = 1, index = 1:sum(pop$n))
  variables$ab_titre$.update()
  variables$inf_time$.update()
  variables$inf_num$.update()

  simulation_loop_safir(
    variables = variables,
    events = events,
    processes = processes,
    timesteps = timesteps,
    variables_dont_update = c("discrete_age", "phase"),
    progress = FALSE
  )

  if (ret_ab) {
    return(variables$ab_titre$get_values())
  } else {
    return(renderer$to_dataframe())
  }

}
