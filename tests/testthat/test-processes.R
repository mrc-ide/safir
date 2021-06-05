test_that("c++ infection process (squire transmission model) returns identical results as R version", {

  # pars
  iso3c <- "GBR"
  pop <- safir:::get_population(iso3c)
  contact_mat <- squire::get_mixing_matrix(iso3c = iso3c)
  pop$n <- as.integer(pop$n / 100)

  parameters <- safir::get_parameters(
    population = pop$n,
    contact_matrix_set = contact_mat,
    iso3c = iso3c,
    time_period = 365
  )
  parameters$beta_set <- parameters$beta_set*rexp(n = length(parameters$beta_set))

  # test states
  n <- 1e5
  dt <- 0.5
  valid_states <- c("S","IMild","ICase","IAsymp")
  state0 <- sample(x = valid_states,size = n,replace = T)
  age0 <- sample.int(n = 17,size = n,replace = T)

  # test R
  states <- individual::CategoricalVariable$new(categories = valid_states,initial_values = state0)
  discrete_age <- individual::IntegerVariable$new(initial_values = age0)
  exposure <- individual::TargetedEvent$new(population_size = n)

  set.seed(5436L)
  inf_proc_R <- infection_process(parameters = parameters,variables = list(states=states,discrete_age=discrete_age),events = list(exposure=exposure),dt =dt)
  inf_proc_R(timestep = 100)

  sched_R <- exposure$get_scheduled()$to_vector()

  # test C++
  states1 <- individual::CategoricalVariable$new(categories = valid_states,initial_values = state0)
  discrete_age1 <- individual::IntegerVariable$new(initial_values = age0)
  exposure1 <- individual::TargetedEvent$new(population_size = n)

  set.seed(5436L)
  inf_proc_Cpp <- infection_process_cpp(parameters = parameters,variables = list(states=states1,discrete_age=discrete_age1),events = list(exposure=exposure1),dt =dt)
  individual:::execute_process(process = inf_proc_Cpp,timestep = 100)

  sched_Cpp <- exposure1$get_scheduled()$to_vector()

  expect_identical(sort(sched_R), sort(sched_Cpp))
})


test_that("c++ infection process (nimue vaccine model) returns identical results as R version", {

  iso3c <- "GBR"
  pop <- safir:::get_population(iso3c)
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
  individual:::execute_process(process = inf_proc_Cpp,timestep = 100)

  sched_Cpp <- exposure1$get_scheduled()$to_vector()

  expect_identical(sort(sched_R), sort(sched_Cpp))
})

# test_that("test create_setup_process", {
#   # create our dummy variables for create_setup_process
#   # these are mock arguments to be passed to the function
#   # Because create_steup_process is istelf a wrapper to a function that takes
#   # an argument api, we will mock the arguments passed to the api
#
#   pop <- get_population("ATG")
#
#   pop$n <- as.integer(pop$n / 100)
#
#   parameters <- get_parameters(
#     iso3c = "ATG",
#     population = pop$n,
#     R0 = 2,
#     time_period = 100,
#     tt_contact_matrix = 0
#   )
#
#   events <- create_events(parameters)
#   variables <- create_variables(pop, parameters)
#   attach_event_listeners(
#     variables,
#     events,
#     parameters
#   )
#
#   events$severe_infection <- mock_event()
#   events$mild_infection <- mock_event()
#   events$asymp_infection <- mock_event()
#
#   variables$states <- mock_category("E", rep("E", 6))
#
#   # Stub delays
#   mockery::stub(create_setup_process, "r_erlang", mockery::mock(c(0.5, 0.5), c(0.4, 0.4), c(0.7, 0.7)))
#   # Stub tree probs
#   mockery::stub(create_setup_process, 'bernoulli_multi_p',
#                 mockery::mock(c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE),
#                               c(TRUE, TRUE, FALSE, FALSE)))
#
#   create_setup_process(
#     parameters,
#     events,
#     variables
#   )
#
#   # First call is to schedule the severe infection for individuals 1 & 2
#   # who schedule in 0.5 + 1 days
#   have_move1 <- individual::Bitset$new(6)
#   have_move1 <- have_move1$insert(1:2)
#   mockery::expect_args(
#     events$severe_infection$schedule,
#     1,
#     have_move1,
#     c(0.5, 0.5)
#   )
#
#   # The second call is to schedule the asymptomatic infection for individuals
#   # 5 & 6 who schedule in 0.4 + 1 days
#   have_move3 <- individual::Bitset$new(6)
#   have_move3 <- have_move3$insert(5:6)
#   mockery::expect_args(
#     events$asymp_infection$schedule,
#     1,
#     have_move3,
#     c(0.4, 0.4)
#   )
#
#   # The thirs call is to schedule the mild infection for individuals 3 & 4
#   # who schedule in 0.7 + 1 days
#   have_move2 <- individual::Bitset$new(6)
#   have_move2 <- have_move2$insert(3:4)
#   mockery::expect_args(
#     events$mild_infection$schedule,
#     1,
#     have_move2,
#     c(0.7, 0.7)
#   )
# })
#
# test_that("test that create_pocesses works", {
#
#   pop <- get_population("AFG")
#   pop$n <- as.integer(pop$n / 10000)
#
#   R0 <- 2
#   time_period <- 1000
#   tt_contact_matrix <- 0
#
#   psq <- get_parameters(
#     iso3c = "AFG",
#     population = pop$n,
#     R0 = R0,
#     time_period = time_period,
#     tt_contact_matrix = tt_contact_matrix
#   )
#
#   variables <- create_variables(pop, psq)
#   events <- create_events(psq)
#   renderer <- individual::Render$new(time_period)
#   # Check create_processes has worked correctly
#   p <- create_processes(  events,
#                           variables,
#                           parameters,
#                           renderer)
#   for (process in p) {
#     expect(is.function(process) || inherits(process, 'externalptr'),
#            'Process is not a function')
#   }
#
# })
