test_that("c++ infection process (multiple doses, no types) returns identical results as R version", {

  library(nimue)

  # pars
  iso3c <- "GBR"
  pop <- safir::get_population(iso3c)
  contact_mat <- squire::get_mixing_matrix(iso3c = iso3c)
  pop$n <- as.integer(pop$n / 100)

  parameters <- safir::get_parameters(
    population = pop$n,
    contact_matrix_set = contact_mat,
    iso3c = iso3c,
    time_period = 365,
    dt = 1
  )
  parameters$beta_set <- parameters$beta_set*rexp(n = length(parameters$beta_set))

  # test states
  n <- 1e5
  dt <- 0.5
  valid_states <- c("S","IMild","ICase","IAsymp")
  state0 <- sample(x = valid_states,size = n,replace = T)
  age0 <- sample.int(n = 17,size = n,replace = T)
  vaccine_parameters <- get_vaccine_ab_titre_parameters(vaccine = "Pfizer")
  ab_titre0 <- log(10^rnorm(n = n, mean = log10(vaccine_parameters$mu_ab[1]),sd = vaccine_parameters$std10))
  parameters <- make_vaccine_parameters(safir_parameters = parameters,vaccine_ab_parameters = vaccine_parameters,vaccine_set = rep(100,365),dose_period = c(NaN, 10),strategy_matrix = nimue::strategy_matrix(strategy = "Elderly"),next_dose_priority_matrix = matrix(0,nrow = 1,ncol = 17))

  # test R
  states <- individual::CategoricalVariable$new(categories = valid_states,initial_values = state0)
  discrete_age <- individual::IntegerVariable$new(initial_values = age0)
  ab_titre <- individual::DoubleVariable$new(initial_values = ab_titre0)
  exposure <- individual::TargetedEvent$new(population_size = n)

  set.seed(5436L)
  inf_proc_R <- infection_process_vaccine(parameters = parameters,variables = list(states=states,discrete_age=discrete_age,ab_titre=ab_titre),events = list(exposure=exposure),dt =dt)
  inf_proc_R(timestep = 100)

  sched_R <- exposure$get_scheduled()$to_vector()

  # test C++
  states1 <- individual::CategoricalVariable$new(categories = valid_states,initial_values = state0)
  discrete_age1 <- individual::IntegerVariable$new(initial_values = age0)
  ab_titre1 <- individual::DoubleVariable$new(initial_values = ab_titre0)
  exposure1 <- individual::TargetedEvent$new(population_size = n)

  set.seed(5436L)
  inf_proc_Cpp <- infection_process_vaccine_cpp(parameters = parameters,variables = list(states=states1,discrete_age=discrete_age1,ab_titre=ab_titre1),events = list(exposure=exposure1),dt =dt)
  execute_process(process = inf_proc_Cpp,timestep = 100)

  sched_Cpp <- exposure1$get_scheduled()$to_vector()

  expect_identical(sort(sched_R), sort(sched_Cpp))
})


test_that("infection process with vaccines: infection_efficacy prevents transmission (S->E) event queueing", {

  library(nimue)
  library(individual)

  # run 1: lots of vaccine

  R0 <- 15
  time_period <- 3
  dt <- 1

  pop <- get_population("GBR")
  pop$n <- rep(1e3, 17) + rep(100, 17)

  parameters <- safir::get_parameters(
    population = pop$n,
    seeding_cases = 0,
    iso3c = "GBR",
    R0 = R0,
    time_period = time_period,
    dt = dt
  )
  parameters$S_0 <- parameters$S_0 - rep(100, 17)
  parameters$ICase1_0 <- rep(100, 17)

  # vaccine parameters
  vaccine_coverage_mat <- nimue::strategy_matrix(strategy = "All",max_coverage = 0.95)
  next_dose_priority <- matrix(nrow = 0, ncol = 17)
  ab_parameters <- get_vaccine_ab_titre_parameters(vaccine = "Moderna", max_dose = 1,correlated = FALSE, max_ab = 10)
  ab_parameters$mu_ab[1] <- 4

  vaccine_doses <- 1
  dose_period <- NaN

  # vaccines
  vaccine_set <- c(1e5, 1e5, 1e5)

  # combine parameters and verify
  parameters <- make_vaccine_parameters(
    safir_parameters = parameters,
    vaccine_ab_parameters = ab_parameters,
    vaccine_set = vaccine_set,
    dose_period = dose_period,
    strategy_matrix = vaccine_coverage_mat,
    next_dose_priority_matrix = next_dose_priority
  )

  # create variables
  timesteps <- parameters$time_period/dt
  variables <- create_variables(pop = pop, parameters = parameters)
  variables <- create_vaccine_variables(variables = variables,parameters = parameters)

  # create events
  events <- create_events(parameters = parameters)
  events <- create_events_vaccination(events = events,parameters = parameters)
  attach_event_listeners(variables = variables,events = events,parameters = parameters, dt = dt)
  attach_event_listeners_vaccination(variables = variables,events = events,parameters = parameters,dt = dt)

  # make renderers
  incidence_renderer <- individual::Render$new(timesteps)

  # processes
  processes <- list(
    vaccine_ab_titre_process(parameters = parameters,variables = variables,dt = dt),
    vaccination_process(parameters = parameters,variables = variables,events = events,dt = dt),
    infection_process_vaccine_cpp(parameters = parameters,variables = variables,events = events,dt = dt)
  )

  attach_tracking_listener_incidence(events = events, renderer = incidence_renderer)
  setup_events(parameters = parameters,events = events,variables = variables,dt = dt)

  # no vaccination prior to sim
  expect_equal(variables$dose_num$get_size_of(set = 0), sum(parameters$population))

  simulation_loop_safir(
    variables = variables,
    events = events,
    processes = processes,
    timesteps = timesteps,
    variables_dont_update = c("discrete_age", "phase"),
    progress = FALSE
  )

  vaccine_run_incidence <- incidence_renderer$to_dataframe()

  expect_true(variables$dose_num$get_size_of(set = 0) < sum(parameters$population))


  # run 2 no vaccines

  # vaccines
  vaccine_set <- vaccine_set*0

  # combine parameters and verify
  parameters <- make_vaccine_parameters(
    safir_parameters = parameters,
    vaccine_ab_parameters = ab_parameters,
    vaccine_set = vaccine_set,
    dose_period = dose_period,
    strategy_matrix = vaccine_coverage_mat,
    next_dose_priority_matrix = next_dose_priority
  )

  # create variables
  timesteps <- parameters$time_period/dt
  variables <- create_variables(pop = pop, parameters = parameters)
  variables <- create_vaccine_variables(variables = variables,parameters = parameters)

  # create events
  events <- create_events(parameters = parameters)
  events <- create_events_vaccination(events = events,parameters = parameters)
  attach_event_listeners(variables = variables,events = events,parameters = parameters, dt = dt)
  attach_event_listeners_vaccination(variables = variables,events = events,parameters = parameters,dt = dt)

  # make renderers
  incidence_renderer <- individual::Render$new(timesteps)

  # processes
  processes <- list(
    vaccine_ab_titre_process(parameters = parameters,variables = variables,dt = dt),
    vaccination_process(parameters = parameters,variables = variables,events = events,dt = dt),
    infection_process_vaccine_cpp(parameters = parameters,variables = variables,events = events,dt = dt)
  )

  attach_tracking_listener_incidence(events = events, renderer = incidence_renderer)
  setup_events(parameters = parameters,events = events,variables = variables,dt = dt)

  # no vaccination prior to sim
  expect_equal(variables$dose_num$get_size_of(set = 0), sum(parameters$population))

  simulation_loop_safir(
    variables = variables,
    events = events,
    processes = processes,
    timesteps = timesteps,
    variables_dont_update = c("discrete_age", "phase"),
    progress = FALSE
  )

  no_vaccine_run_incidence <- incidence_renderer$to_dataframe()

  expect_equal(variables$dose_num$get_size_of(set = 0), sum(parameters$population))

  # tests
  expect_true(all(no_vaccine_run_incidence[-1, "incidence"] > vaccine_run_incidence[-1, "incidence"]))

})
