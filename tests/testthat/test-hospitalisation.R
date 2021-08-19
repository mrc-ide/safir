test_that('allocate_treatment can allocate when limit is not exceeded', {
  pop <- get_population("ATG")

  pop$n <- as.integer(pop$n / 100)

  parameters <- get_parameters(
    iso3c = "ATG",
    population = pop$n,
    R0 = 2,
    time_period = 100,
    tt_contact_matrix = 0,
    dt = 1
  )

  events <- create_events(parameters)
  variables <- create_variables(pop, parameters)
  attach_event_listeners(
    variables,
    events,
    parameters,
    dt = 0.1
  )

  variables$states <- mock_category(c("S", 'IOxGetDie', 'IOxGetLive', 'IRec'),
                                    c("S", "S", "S", "S", "S"))

  to_treat <- individual::Bitset$new(5)
  to_treat <- to_treat$insert(1:5)
  allocated <- allocate_treatment(
    variables,
    to_treat,
    c('IOxGetDie', 'IOxGetLive', 'IRec'),
    100
  )

  expect_equal(allocated$size(), 5)
})

test_that('allocate_treatment can allocate when limit is exceeded', {
  pop <- get_population("ATG")

  pop$n <- as.integer(pop$n / 100)

  parameters <- get_parameters(
    iso3c = "ATG",
    population = pop$n,
    R0 = 2,
    time_period = 100,
    tt_contact_matrix = 0,
    dt = 1
  )

  events <- create_events(parameters)
  variables <- create_variables(pop, parameters)
  attach_event_listeners(
    variables,
    events,
    parameters,
    dt = 0.1
  )

  variables$states <- mock_category(c("S", 'IOxGetDie', 'IOxGetLive', 'IRec'),
                                    c("S", "S", "IRec", "IRec", "IRec"))

  to_treat <- individual::Bitset$new(5)
  to_treat <- to_treat$insert(1:2)
  allocated <- allocate_treatment(
    variables,
    to_treat,
    c('IOxGetDie', 'IOxGetLive', 'IRec'),
    3
  )

  expect_equal(allocated$size(), 0)
})
