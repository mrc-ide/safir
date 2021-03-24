test_that("create_event_based_processes assigns a listener to each event", {

  pop <- get_population("ATG")

  pop$n <- as.integer(pop$n / 100)

  psq <- get_parameters(
    iso3c = "ATG",
    population = pop$n,
    R0 = 2,
    time_period = 100,
    tt_contact_matrix = 0
  )

  events <- create_events()
  states <- create_states(psq)
  variables <- create_variables(pop, psq)
  human <- create_human(states, variables, events)

  create_event_based_processes(
    human,
    states,
    variables,
    events,
    psq
  )

  for (event in events) {
    expect_gt(length(event$listeners), 0)
  }
})

test_that("test create_infection_update_listener", {

  human <- mockery::mock()
  IMild <- mockery::mock()

  ret <- create_infection_update_listener(
    human,
    IMild
  )

  api <- list(queue_state_update = mockery::mock())
  to_move <- mockery::mock()
  ret(api, to_move)
  mockery::expect_args(api$queue_state_update, 1, human, IMild, to_move)

})

test_that("test create_progression_listener", {

  events <- list(exposure = mockery::mock())
  shift <- 0
  duration <- 1
  target <- mockery::mock()
  func_mock <- mockery::mock()
  r_erlang_mock <- mockery::mock(c(TRUE, TRUE, TRUE, TRUE))

  ret <- create_progression_listener(
    events$exposure,
    duration,
    shift,
    func_mock
  )

  api <- list(schedule = mockery::mock())

  mockery::stub(ret, 'r_erlang', mockery::mock(c(TRUE, TRUE, TRUE, TRUE)))

  with_mock(
    'safir::r_erlang' = r_erlang_mock,
    ret2 <- api$schedule(
      events$exposure,
      target,
      rep(.5, 4)
    )
  )

  ret(api, target)

  mockery::expect_args(
    api$schedule,
    1,
    event = events$exposure,
    target = target,
    func_mock = c(0.5, 0.5, 0.5, 0.5)
  )

  mockery::expect_called(api$schedule, 2)

})

test_that("test create_exposure_update_listener", {

  human <- mockery::mock()
  states <- mockery::mock()

  ret <- create_exposure_update_listener(
    human,
    states,
    events,
    variables,
    parameters
  )

  api <- list(get_variable = mockery::mock(), schedule = mockery::mock())

  variables <- list(discrete_age = mockery::mock())

  events <- list(exposure = mockery::mock(),
                 mild_infection = mockery::mock(),
                 asymp_infection = mockery::mock(),
                 severe_infection = mockery::mock())

  problist <- seq(1,17,1)
  exposure = mockery::mock()

  parameters <- list(dur_E = mockery::mock(),
                     prob_hosp = list(mockery::mock(problist)),
                     prob_asymp = list(mockery::mock(problist)))

  to_move <- c(6, 7, 3, 2, 5, 8)

  mockery::stub(ret, 'r_erlang', mockery::mock(c(0.5, 0.5), c(0.4, 0.4),
                                               c(0.7, 0.7)))

  mockery::stub(ret, 'bernoulli_multi_p',
                mockery::mock(c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE),
                              c(TRUE, TRUE, FALSE, FALSE)))
  ret(api, to_move)

  # The process should have called api$schedule three times

  # First call is to schedule the severe infection for individuals 6, 7
  # who schedule in 0 days
  mockery::expect_args(
    api$schedule,
    1,
    event = events$severe_infection,
    target = c(6, 7),
    delay = c(1.5, 1.5)
  )

  # The second call is to schedule the mild infection for individuals 3, 2
  # who schedule in 0 days
  mockery::expect_args(
    api$schedule,
    2,
    event = events$mild_infection,
    target = c(3, 2),
    delay = c(1.4, 1.4)
  )

  # The third call is to schedule the asymptomatic infection for individuals
  # 5, 8 who schedule in 0 days
  mockery::expect_args(
    api$schedule,
    3,
    event = events$asymp_infection,
    target = c(5, 8),
    delay = c(1.7, 1.7)
  )

})


test_that("test initialise_progression", {

  event <- mockery::mock()
  human <- mockery::mock()
  from_state <- mockery::mock()
  duration <- 1
  r_erlang_mock <- mockery::mock(c(TRUE, TRUE, TRUE, TRUE))

  ret <- initialise_progression(event, human, from_state, duration)

  api <- list(schedule = mockery::mock(), get_state = mockery::mock())
  target <- mockery::mock()

  ret(api, target)

  mockery::expect_args(api$get_state, 1, human, from_state)

  with_mock(
    'safir::r_erlang' = r_erlang_mock,
    ret3 <- api$schedule(
      event,
      target,
      rep(.5, 4)
    )
  )

})
