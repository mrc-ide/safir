test_that("c++ infection process (nimue vaccine model) returns identical results as R version", {

  iso3c <- "GBR"
  pop <- safir::get_population(iso3c)
  pop$n <- as.integer(pop$n / 2e2)
  contact_mat <- squire::get_mixing_matrix(iso3c = iso3c)

  tmax <- 250
  R0 <- 4
  dt <- 0.01

  parameters <- get_parameters_nimue(
    population = pop$n,
    contact_mat = contact_mat,
    time_period = tmax,
    R0 = R0,
    max_vaccine = c(0, seq(10, 1e4, length.out = 40)),
    tt_vaccine = c(0, seq(10, 50, length.out = 40)),
    vaccine_efficacy_disease = rep(0, 17),
    vaccine_efficacy_infection = rep(0.9, 17)
  )

  # test states
  n <- 1e5
  dt <- 0.5
  valid_states <- c("S","IMild","ICase","IAsymp")
  state0 <- sample(x = valid_states,size = n,replace = T)
  age0 <- sample.int(n = 17,size = n,replace = T)
  vaxx0 <- sample.int(n = 4,size = n,replace = T)

  # test R
  states <- individual::CategoricalVariable$new(categories = valid_states,initial_values = state0)
  discrete_age <- individual::IntegerVariable$new(initial_values = age0)
  vaccine_states <- individual::IntegerVariable$new(initial_values = vaxx0)
  exposure <- individual::TargetedEvent$new(population_size = n)

  set.seed(5436L)
  inf_proc_R <- infection_process_nimue(parameters = parameters,variables = list(states=states,vaccine_states=vaccine_states,discrete_age=discrete_age),events = list(exposure=exposure),dt =dt)
  inf_proc_R(timestep = 100)

  sched_R <- exposure$get_scheduled()$to_vector()

  # test C++
  states1 <- individual::CategoricalVariable$new(categories = valid_states,initial_values = state0)
  discrete_age1 <- individual::IntegerVariable$new(initial_values = age0)
  vaccine_states1 <- individual::IntegerVariable$new(initial_values = vaxx0)
  exposure1 <- individual::TargetedEvent$new(population_size = n)

  set.seed(5436L)
  inf_proc_Cpp <- infection_process_nimue_cpp(parameters = parameters,variables = list(states=states1,vaccine_states=vaccine_states1,discrete_age=discrete_age1),events = list(exposure=exposure1),dt =dt)
  execute_process(process = inf_proc_Cpp,timestep = 100)

  sched_Cpp <- exposure1$get_scheduled()$to_vector()

  expect_identical(sort(sched_R), sort(sched_Cpp))
})

test_that("c++ vaccination process (nimue vaccine model) returns identical results as R version", {

  iso3c <- "GBR"
  pop <- safir::get_population(iso3c)
  pop$n <- as.integer(pop$n / 2e2)
  contact_mat <- squire::get_mixing_matrix(iso3c = iso3c)

  tmax <- 250
  R0 <- 4
  dt <- 0.01
  tt_vaccine <- c(0, 10:150)
  max_vaccine <- c(0, seq(1e3, 5e4, length.out = length(tt_vaccine)-1))

  parameters <- get_parameters_nimue(
    population = pop$n,
    contact_mat = contact_mat,
    time_period = tmax,
    tt_vaccine=tt_vaccine,
    max_vaccine=max_vaccine,
    R0 = R0,
  )

  # test states
  n <- 1e5
  dt <- 0.5
  valid_states <- c("S","E","R","IMild","ICase","IAsymp")
  state0 <- sample(x = valid_states,size = n,replace = T)
  age0 <- sample.int(n = 17,size = n,replace = T)
  vaxx0 <- sample.int(n = n,size = 2e4,replace = T)

  # test R
  states <- individual::CategoricalVariable$new(categories = valid_states,initial_values = state0)
  discrete_age <- individual::IntegerVariable$new(initial_values = age0)
  eligible <- individual::Bitset$new(n)
  empty <- individual::Bitset$new(n)
  vaccinated <- individual::Bitset$new(n)$insert(vaxx0)
  v0_to_v1v2 <- individual::TargetedEvent$new(population_size = n)

  vars <- list(
    states=states,
    discrete_age=discrete_age,
    eligible=eligible,
    empty=empty,
    vaccinated=vaccinated
  )

  set.seed(959256L)
  inf_proc_R <- vaccination_process_nimue(
    parameters = parameters,
    variables = vars,
    events = list(v0_to_v1v2=v0_to_v1v2),
    dt =dt
  )
  inf_proc_R(timestep = 120)

  sched_R <- v0_to_v1v2$get_scheduled()$to_vector()

  # test C++
  states1 <- individual::CategoricalVariable$new(categories = valid_states,initial_values = state0)
  discrete_age1 <- individual::IntegerVariable$new(initial_values = age0)
  eligible1 <- individual::Bitset$new(n)
  empty1 <- individual::Bitset$new(n)
  vaccinated1 <- individual::Bitset$new(n)$insert(vaxx0)
  v0_to_v1v21 <- individual::TargetedEvent$new(population_size = n)

  vars <- list(
    states=states1,
    discrete_age=discrete_age1,
    eligible=eligible1,
    empty=empty1,
    vaccinated=vaccinated1
  )

  set.seed(959256L)
  inf_proc_Cpp <- vaccination_process_nimue_cpp(
    parameters = parameters,
    variables = vars,
    events = list(v0_to_v1v2=v0_to_v1v21),
    dt =dt
  )
  execute_process(process = inf_proc_Cpp,timestep = 120)

  sched_Cpp <- v0_to_v1v21$get_scheduled()$to_vector()

  expect_identical(sort(sched_R), sort(sched_Cpp))
})
