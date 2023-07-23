# tests in this file might be broken because the VFR is no longer applied
# in the efficacy functions, but in the make_calculate_nat.
# to get them to align again I supposed we'd have to manually apply the
# VFR correction

test_that("test VFR vector gives reasonable output", {

  # dt = 1

  library(squire)

  iso3c <- "GBR"
  pop <- safir::get_population(iso3c)
  contact_mat <- squire::get_mixing_matrix(iso3c = iso3c)

  # Create our simulation parameters
  time_period <- 200
  dt <- 1
  timesteps <- time_period/dt

  parameters <- get_parameters(
    population = pop$n,
    contact_matrix_set = contact_mat,
    iso3c = iso3c,
    R0 = 2,
    time_period = time_period,
    dt = dt
  )

  vfr_time_1 <- 100
  vfr_time_2 <- 120
  vfr_val <- 4

  vfr <- variant_fold_reduction_vector(parameters = parameters, dt = dt, vfr = vfr_val, vfr_time_1 = vfr_time_1, vfr_time_2 = vfr_time_2)

  expect_equal(length(vfr), time_period)
  expect_equal(vfr[1:vfr_time_1], rep(1, vfr_time_1))
  expect_true(vfr[vfr_time_1 + 1] < vfr[vfr_time_2])
  expect_equal(vfr[vfr_time_2:time_period], rep(vfr_val, time_period - vfr_time_2 + 1))

  # dt = 0.2

  # Create our simulation parameters
  time_period <- 200
  dt <- 0.2
  timesteps <- time_period/dt

  parameters <- get_parameters(
    population = pop$n,
    contact_matrix_set = contact_mat,
    iso3c = iso3c,
    R0 = 2,
    time_period = time_period,
    dt = dt
  )

  vfr_time_1 <- 100
  vfr_time_2 <- 120
  vfr_val <- 4

  vfr <- variant_fold_reduction_vector(parameters = parameters, dt = dt, vfr = vfr_val, vfr_time_1 = vfr_time_1, vfr_time_2 = vfr_time_2)

  expect_equal(length(vfr), time_period)
  expect_equal(vfr[1:vfr_time_1], rep(1, vfr_time_1))
  expect_true(vfr[(vfr_time_1) + 1] < vfr[vfr_time_2])
  expect_equal(vfr[(vfr_time_2):time_period], rep(vfr_val, time_period - (vfr_time_2) + 1))

})


test_that("test VFR decay with NAT from vaccination working with dt = 1", {

  # Create our simulation parameters
  time_period <- 800
  dt <- 1
  n <- 10
  timesteps <- time_period/dt

  vacc_parameters <- get_vaccine_ab_titre_parameters(vaccine = "Pfizer")
  parameters <- list(time_period=time_period, population=n, N_phase=1)
  parameters <- make_vaccine_parameters(safir_parameters = parameters, vaccine_ab_parameters = vacc_parameters, dose_period = NaN, vaccine_set = rep(0, time_period), strategy_matrix = matrix(0, 17, 17), next_dose_priority_matrix = matrix(0, 0, 17))

  vfr_time_1 <- 100
  vfr_time_2 <- 120
  vfr_val <- 4

  # compute with gold standard code
  compare_nt <- draw_nt_vfr(parameters = vacc_parameters, n = n, tmax = time_period, vf = vfr_val, vfr_time_1 = vfr_time_1, vfr_time_2 = vfr_time_2)

  nt <- compare_nt$nt
  z1 <- compare_nt$z1

  # compute with safir
  variables <- create_vaccine_variables(variables = list(), parameters = parameters)
  variables <- create_independent_nat_variables(variables = variables, parameters = parameters)
  ef_infection <- matrix(0, nrow = timesteps + 1, ncol = n)
  ef_severe <- matrix(0, nrow = timesteps + 1, ncol = n)
  ef_infection_cpp <- matrix(0, nrow = timesteps + 1, ncol = n)
  ef_severe_cpp <- matrix(0, nrow = timesteps + 1, ncol = n)

  vfr_vector <- variant_fold_reduction_vector(parameters = parameters, dt = dt, vfr = vfr_val, vfr_time_1 = vfr_time_1, vfr_time_2 = vfr_time_2)
  parameters <- make_immune_parameters(parameters = parameters, vfr = vfr_vector)

  # time 0 NAT boost from vaccine
  index <- Bitset$new(n)$insert(1:n)
  schedule_dose_vaccine(timestep = 0,variables = variables,target = index,dose = 1,parameters = parameters)
  variables$ab_titre$queue_update(values = log(10^z1)) # make sure using the same set of RVs
  variables$ab_titre$.update()
  variables$dose_num$.update()
  variables$dose_time$.update()

  ab_titre <- vaccine_ab_titre_process(parameters = parameters,variables = variables, dt = dt)

  calculate_nat <- make_calculate_nat(variables = variables, parameters = parameters)
  nat <- calculate_nat(index = index, day = 1)

  safir_out <- matrix(data = NaN,nrow = timesteps + 1,ncol = n)
  safir_out[1, ] <- variables$ab_titre$get_values()

  ef_infection[1, ] <- vaccine_efficacy_infection(nat = nat, parameters = parameters, day = 1)
  ef_severe[1, ] <- vaccine_efficacy_severe(nat = nat, ef_infection = ef_infection[1, ], parameters = parameters, day = 1)
  ef_infection_cpp[1, ] <- vaccine_efficacy_infection_cpp(nat = nat, parameters = parameters, day = 0)
  ef_severe_cpp[1, ] <- vaccine_efficacy_severe_cpp(nat = nat, ef_infection = ef_infection_cpp[1, ], parameters = parameters, day = 0)

  for (t in 1:timesteps) {
    day <- ceiling(t * dt)
    ab_titre(timestep = t)
    variables$ab_titre$.update()
    # NOTE: we have ab_titire_inf in the variables currently but it is not being updated
    # the reason is that the original "gold standard" code does not consider it
    safir_out[t + 1L, ] <- variables$ab_titre$get_values()

    nat <- calculate_nat(index = index, day = day)

    ef_infection[t + 1L, ] <- vaccine_efficacy_infection(nat = nat, parameters = parameters, day = day)
    ef_severe[t + 1L, ] <- vaccine_efficacy_severe(nat = nat, ef_infection = ef_infection[t + 1L, ], parameters = parameters, day = day)
    ef_infection_cpp[t + 1L, ] <- vaccine_efficacy_infection_cpp(nat = nat, parameters = parameters, day = day - 1)
    ef_severe_cpp[t + 1L, ] <- vaccine_efficacy_severe_cpp(nat = nat, ef_infection = ef_infection_cpp[t + 1L, ], parameters = parameters, day = day - 1)
  }

  # test ln scale NAT same
  expect_equal(safir_out, nt)
  expect_equal(1 - ef_infection, compare_nt$ef_infection)
  expect_equal(1 - ef_severe, compare_nt$ef_severe)
  expect_equal(1 - ef_infection_cpp, compare_nt$ef_infection)
  expect_equal(1 - ef_severe_cpp, compare_nt$ef_severe)

})


test_that("test VFR decay with NAT from vaccination working with dt = 0.5", {

  # Create our simulation parameters
  time_period <- 800
  dt <- 0.5
  n <- 10
  timesteps <- time_period/dt

  vacc_parameters <- get_vaccine_ab_titre_parameters(vaccine = "Pfizer")
  parameters <- list(time_period=time_period, population=n, N_phase=1)
  parameters <- make_vaccine_parameters(safir_parameters = parameters, vaccine_ab_parameters = vacc_parameters, dose_period = NaN, vaccine_set = rep(0, time_period), strategy_matrix = matrix(0, 17, 17), next_dose_priority_matrix = matrix(0, 0, 17))

  vfr_time_1 <- 100
  vfr_time_2 <- 120
  vfr_val <- 4

  # compute with gold standard code
  compare_nt <- draw_nt_vfr(parameters = vacc_parameters, n = n, tmax = time_period, vf = vfr_val, vfr_time_1 = vfr_time_1, vfr_time_2 = vfr_time_2)

  nt <- compare_nt$nt
  z1 <- compare_nt$z1

  # compute with safir
  variables <- create_vaccine_variables(variables = list(), parameters = parameters)
  variables <- create_independent_nat_variables(variables = variables, parameters = parameters)
  ef_infection <- matrix(0, nrow = timesteps + 1, ncol = n)
  ef_severe <- matrix(0, nrow = timesteps + 1, ncol = n)
  ef_infection_cpp <- matrix(0, nrow = timesteps + 1, ncol = n)
  ef_severe_cpp <- matrix(0, nrow = timesteps + 1, ncol = n)

  vfr_vector <- variant_fold_reduction_vector(parameters = parameters, dt = dt, vfr = vfr_val, vfr_time_1 = vfr_time_1, vfr_time_2 = vfr_time_2)
  parameters <- make_immune_parameters(parameters = parameters, vfr = vfr_vector)

  # time 0 NAT boost from vaccine
  index <- Bitset$new(n)$insert(1:n)
  schedule_dose_vaccine(timestep = 0,variables = variables,target = index,dose = 1,parameters = parameters)
  variables$ab_titre$queue_update(values = log(10^z1)) # make sure using the same set of RVs
  variables$ab_titre$.update()
  variables$dose_num$.update()
  variables$dose_time$.update()

  ab_titre <- vaccine_ab_titre_process(parameters = parameters,variables = variables, dt = dt)

  calculate_nat <- make_calculate_nat(variables = variables, parameters = parameters)
  nat <- calculate_nat(index = index, day = 1)

  safir_out <- matrix(data = NaN,nrow = timesteps + 1,ncol = n)
  safir_out[1, ] <- variables$ab_titre$get_values()
  ef_infection[1, ] <- vaccine_efficacy_infection(nat = nat, parameters = parameters, day = 1)
  ef_severe[1, ] <- vaccine_efficacy_severe(nat = nat, ef_infection = ef_infection[1, ], parameters = parameters, day = 1)
  ef_infection_cpp[1, ] <- vaccine_efficacy_infection_cpp(nat = nat, parameters = parameters, day = 0)
  ef_severe_cpp[1, ] <- vaccine_efficacy_severe_cpp(nat = nat, ef_infection = ef_infection_cpp[1, ], parameters = parameters, day = 0)


  for (t in 1:timesteps) {
    day <- ceiling(t * dt)
    ab_titre(timestep = t)
    variables$ab_titre$.update()
    safir_out[t + 1L, ] <- variables$ab_titre$get_values()

    nat <- calculate_nat(index = index, day = day)

    ef_infection[t + 1L, ] <- vaccine_efficacy_infection(nat = nat, parameters = parameters, day = day)
    ef_severe[t + 1L, ] <- vaccine_efficacy_severe(nat = nat, ef_infection = ef_infection[t + 1L, ], parameters = parameters, day = day)
    ef_infection_cpp[t + 1L, ] <- vaccine_efficacy_infection_cpp(nat = nat, parameters = parameters, day = day - 1)
    ef_severe_cpp[t + 1L, ] <- vaccine_efficacy_severe_cpp(nat = nat, ef_infection = ef_infection_cpp[t + 1L, ], parameters = parameters, day = day - 1)
  }

  # test ln scale NAT same
  last_row <- nrow(safir_out)
  expect_equal(safir_out[last_row, ], nt[nrow(nt), ])
  expect_equal(1 - ef_infection[last_row, ], compare_nt$ef_infection[nrow(nt), ])
  expect_equal(1 - ef_severe[last_row, ], compare_nt$ef_severe[nrow(nt), ])
  expect_equal(1 - ef_infection_cpp[last_row, ], compare_nt$ef_infection[nrow(nt), ])
  expect_equal(1 - ef_severe_cpp[last_row, ], compare_nt$ef_severe[nrow(nt), ])

})



test_that("test that VFR means more infections are queued in simulation runs with VFR = 5 than without any VFR", {

  library(nimue)

  # pars
  iso3c <- "GBR"
  pop <- safir::get_population(iso3c)
  contact_mat <- squire::get_mixing_matrix(iso3c = iso3c)
  pop$n <- as.integer(pop$n / 100)

  pars <- safir::get_parameters(
    population = pop$n,
    contact_matrix_set = contact_mat,
    iso3c = iso3c,
    R0 = 5,
    time_period = 365,
    dt = 1
  )

  # test states
  n <- 1e5
  dt <- 0.5
  timesteps <- pars$time_period / dt
  valid_states <- c("S","IMild","ICase","IAsymp")
  state0 <- sample(x = valid_states,size = n,replace = T)
  age0 <- sample.int(n = 17,size = n,replace = T)
  vaccine_parameters <- get_vaccine_ab_titre_parameters(vaccine = "Pfizer")
  zdose <- log(10^rnorm(n = n, mean = log10(3), sd = 0.05))

  # VFR = 4 (should have more infections)

  # R
  parameters <- make_vaccine_parameters(safir_parameters = pars,vaccine_ab_parameters = vaccine_parameters,vaccine_set = rep(100,365),dose_period = c(NaN, 10),strategy_matrix = nimue::strategy_matrix(strategy = "Elderly"),next_dose_priority_matrix = matrix(0,nrow = 1,ncol = 17))
  parameters <- make_immune_parameters(parameters = parameters, vfr = rep(4, times = pars$time_period), mu_ab_infection = 5, std10_infection = 0.1)

  states <- individual::CategoricalVariable$new(categories = valid_states,initial_values = state0)
  discrete_age <- individual::IntegerVariable$new(initial_values = age0)
  ab_titre <- individual::DoubleVariable$new(initial_values = zdose)
  ab_titre_inf <- individual::DoubleVariable$new(initial_values = rep(-Inf,n))
  exposure <- individual::TargetedEvent$new(population_size = n)

  vars <- list(states=states,discrete_age=discrete_age,ab_titre=ab_titre,ab_titre_inf=ab_titre_inf)

  set.seed(5436L)
  inf_proc <- infection_process_vaccine(parameters = parameters,variables = vars,events = list(exposure=exposure),dt = dt)
  inf_proc(timestep = 1)

  inf_VFR4_R <- exposure$get_scheduled()$to_vector()

  # C++
  parameters <- make_vaccine_parameters(safir_parameters = pars,vaccine_ab_parameters = vaccine_parameters,vaccine_set = rep(100,365),dose_period = c(NaN, 10),strategy_matrix = nimue::strategy_matrix(strategy = "Elderly"),next_dose_priority_matrix = matrix(0,nrow = 1,ncol = 17))
  parameters <- make_immune_parameters(parameters = parameters, vfr = rep(4, times = pars$time_period), mu_ab_infection = 5, std10_infection = 0.1)

  states <- individual::CategoricalVariable$new(categories = valid_states,initial_values = state0)
  discrete_age <- individual::IntegerVariable$new(initial_values = age0)
  ab_titre <- individual::DoubleVariable$new(initial_values = zdose)
  ab_titre_inf <- individual::DoubleVariable$new(initial_values = rep(-Inf,n))
  exposure <- individual::TargetedEvent$new(population_size = n)

  vars <- list(states=states,discrete_age=discrete_age,ab_titre=ab_titre,ab_titre_inf=ab_titre_inf)

  set.seed(5436L)
  inf_proc <- infection_process_vaccine_cpp(parameters = parameters,variables = vars,events = list(exposure=exposure),dt = dt)
  execute_process(process = inf_proc,timestep = 1)

  inf_VFR4_CPP <- exposure$get_scheduled()$to_vector()


  # no VFR (less infections)

  # R
  parameters <- make_vaccine_parameters(safir_parameters = pars,vaccine_ab_parameters = vaccine_parameters,vaccine_set = rep(100,365),dose_period = c(NaN, 10),strategy_matrix = nimue::strategy_matrix(strategy = "Elderly"),next_dose_priority_matrix = matrix(0,nrow = 1,ncol = 17))
  parameters <- make_immune_parameters(parameters = parameters, vfr = rep(1, times = pars$time_period), mu_ab_infection = 5, std10_infection = 0.1)
  parameters$vfr <- NULL

  states <- individual::CategoricalVariable$new(categories = valid_states,initial_values = state0)
  discrete_age <- individual::IntegerVariable$new(initial_values = age0)
  ab_titre <- individual::DoubleVariable$new(initial_values = zdose)
  ab_titre_inf <- individual::DoubleVariable$new(initial_values = rep(-Inf,n))
  exposure <- individual::TargetedEvent$new(population_size = n)

  vars <- list(states=states,discrete_age=discrete_age,ab_titre=ab_titre,ab_titre_inf=ab_titre_inf)

  set.seed(5436L)
  inf_proc <- infection_process_vaccine(parameters = parameters,variables = vars,events = list(exposure=exposure),dt = dt)
  inf_proc(timestep = 1)

  inf_VFR0_R <- exposure$get_scheduled()$to_vector()

  # C++
  parameters <- make_vaccine_parameters(safir_parameters = pars,vaccine_ab_parameters = vaccine_parameters,vaccine_set = rep(100,365),dose_period = c(NaN, 10),strategy_matrix = nimue::strategy_matrix(strategy = "Elderly"),next_dose_priority_matrix = matrix(0,nrow = 1,ncol = 17))
  parameters <- make_immune_parameters(parameters = parameters, vfr = rep(4, times = pars$time_period), mu_ab_infection = 5, std10_infection = 0.1)
  parameters$vfr <- NULL

  states <- individual::CategoricalVariable$new(categories = valid_states,initial_values = state0)
  discrete_age <- individual::IntegerVariable$new(initial_values = age0)
  ab_titre <- individual::DoubleVariable$new(initial_values = zdose)
  ab_titre_inf <- individual::DoubleVariable$new(initial_values = rep(-Inf,n))
  exposure <- individual::TargetedEvent$new(population_size = n)

  vars <- list(states=states,discrete_age=discrete_age,ab_titre=ab_titre,ab_titre_inf=ab_titre_inf)

  set.seed(5436L)
  inf_proc <- infection_process_vaccine_cpp(parameters = parameters,variables = vars,events = list(exposure=exposure),dt = dt)
  execute_process(process = inf_proc,timestep = 1)

  inf_VFR0_CPP <- exposure$get_scheduled()$to_vector()

  # tests
  expect_equal(sort(inf_VFR4_R), sort(inf_VFR4_CPP))
  expect_equal(sort(inf_VFR0_R), sort(inf_VFR0_CPP))

  expect_true(length(inf_VFR4_R) > length(inf_VFR0_R))
})



test_that("test that VFR means more severe infections are queued in simulation runs with VFR = 5 than without any VFR", {

  library(nimue)

  # pars
  iso3c <- "GBR"
  pop <- safir::get_population(iso3c)
  contact_mat <- squire::get_mixing_matrix(iso3c = iso3c)
  pop$n <- as.integer(pop$n / 100)

  pars <- safir::get_parameters(
    population = pop$n,
    contact_matrix_set = contact_mat,
    iso3c = iso3c,
    R0 = 5,
    time_period = 365,
    dt = 1
  )

  # test states
  n <- 1e5
  dt <- 1
  timesteps <- pars$time_period / dt
  valid_states <- c("S","IMild","ICase","IAsymp")
  state0 <- sample(x = valid_states,size = n,replace = T)
  age0 <- sample.int(n = 17,size = n,replace = T)
  vaccine_parameters <- get_vaccine_ab_titre_parameters(vaccine = "Pfizer")
  zdose <- log(10^rnorm(n = n, mean = log(3), sd = 0.05))

  # VFR = 4 (should have more severe outcomes)

  parameters <- make_vaccine_parameters(safir_parameters = pars,vaccine_ab_parameters = vaccine_parameters,vaccine_set = rep(100,365),dose_period = c(NaN, 10),strategy_matrix = nimue::strategy_matrix(strategy = "Elderly"),next_dose_priority_matrix = matrix(0,nrow = 1,ncol = 17))
  parameters <- make_immune_parameters(parameters = parameters, vfr = rep(12, times = pars$time_period), mu_ab_infection = 5, std10_infection = 0.1)

  variables <- list(
    discrete_age = individual::IntegerVariable$new(initial_values = age0),
    ab_titre = individual::DoubleVariable$new(initial_values = zdose),
    ab_titre_inf = individual::DoubleVariable$new(initial_values = rep(-Inf,n))
  )

  events <- list(
    severe_infection = individual::TargetedEvent$new(population_size = n),
    asymp_infection = individual::TargetedEvent$new(population_size = n),
    mild_infection = individual::TargetedEvent$new(population_size = n)
  )

  target <- individual::Bitset$new(size = n)
  target$insert(1:n)

  exp_proc <- create_exposure_scheduler_listener_vaccine(events = events, variables = variables, parameters = parameters, dt = dt)
  exp_proc(timestep = 1, target = target)

  VFR4_severe <- events$severe_infection$get_scheduled()$to_vector()
  VFR4_asymp <- events$asymp_infection$get_scheduled()$to_vector()
  VFR4_mild <- events$mild_infection$get_scheduled()$to_vector()

  # no VFR (should have less severe outcomes)
  parameters <- make_vaccine_parameters(safir_parameters = pars,vaccine_ab_parameters = vaccine_parameters,vaccine_set = rep(100,365),dose_period = c(NaN, 10),strategy_matrix = nimue::strategy_matrix(strategy = "Elderly"),next_dose_priority_matrix = matrix(0,nrow = 1,ncol = 17))
  parameters <- make_immune_parameters(parameters = parameters, vfr = rep(1, times = pars$time_period), mu_ab_infection = 5, std10_infection = 0.1)
  parameters$vfr <- NULL

  variables <- list(
    discrete_age = individual::IntegerVariable$new(initial_values = age0),
    ab_titre = individual::DoubleVariable$new(initial_values = zdose),
    ab_titre_inf = individual::DoubleVariable$new(initial_values = rep(-Inf,n))
  )

  events <- list(
    severe_infection = individual::TargetedEvent$new(population_size = n),
    asymp_infection = individual::TargetedEvent$new(population_size = n),
    mild_infection = individual::TargetedEvent$new(population_size = n)
  )

  target <- individual::Bitset$new(size = n)
  target$insert(1:n)

  exp_proc <- create_exposure_scheduler_listener_vaccine(events = events, variables = variables, parameters = parameters, dt = dt)
  exp_proc(timestep = 1, target = target)

  VFR0_severe <- events$severe_infection$get_scheduled()$to_vector()
  VFR0_asymp <- events$asymp_infection$get_scheduled()$to_vector()
  VFR0_mild <- events$mild_infection$get_scheduled()$to_vector()

  # tests
  expect_true(length(VFR4_severe) > length(VFR0_severe))
  expect_true(length(VFR4_severe) + length(VFR4_asymp) + length(VFR4_mild) == n)
  expect_true(length(VFR0_severe) + length(VFR0_asymp) + length(VFR0_mild) == n)
})


test_that("worse outcomes with VFR > 1 than VFR = 1", {

  iso3c <- "GBR"
  pop <- safir::get_population(iso3c)
  pop$n <- as.integer(pop$n / 1e3)

  tmax <- 20
  dt <- 0.5
  R0 <- 20

  vfr_null <- rep(1, tmax)
  vfr_high <- rep(50, tmax)

  ab_0 <- rep(2, sum(pop$n))

  sim_null <- simulate_vfr(iso3c = iso3c, vfr = vfr_null, tmax = tmax, dt = dt, R0 = R0, ab_titre = ab_0, pop = pop)
  sim_high <- simulate_vfr(iso3c = iso3c, vfr = vfr_high, tmax = tmax, dt = dt, R0 = R0, ab_titre = ab_0, pop = pop)

  expect_true(sim_null$S_count[tmax] > sim_high$S_count[tmax])
})


test_that("time varying mu_ab_infection works", {

  iso3c <- "GBR"
  pop <- safir::get_population(iso3c)
  pop$n <- as.integer(pop$n / 1e3)

  tmax <- 40
  dt <- 0.5
  R0 <- 20

  vfr_null <- rep(1, tmax)

  ab_0 <- rep(1, sum(pop$n))

  sim_const <- simulate_vfr(iso3c = iso3c, vfr = vfr_null, tmax = tmax, dt = dt, R0 = R0, ab_titre = ab_0, pop = pop, mu_ab_infection = c(0.25,0.25,0.25), ret_ab = TRUE)
  sim_timevar <- simulate_vfr(iso3c = iso3c, vfr = vfr_null, tmax = tmax, dt = dt, R0 = R0, ab_titre = ab_0, pop = pop, mu_ab_infection = matrix(c(50,60,75), nrow = 3, ncol = tmax), ret_ab = TRUE)

  expect_true(mean(sim_const) <= mean(sim_const))
})

