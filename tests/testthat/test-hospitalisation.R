test_that('prob_outcome draws correct probabilities', {
  expect_equal(
    prob_outcome(
      target = c(5, 6, 7),
      age = c(1, 1, 1, 2, 3, 4, 5),
      probs = c(.1, .2, .3, .4, .5, .6, .7)
    ),
    c(.3, .4, .5)
  )
})

test_that('schedule_outcome can schedule both successes and failures', {
  api <- list(schedule = mockery::mock())
  success_event <- mockery::mock()
  failure_event <- mockery::mock()
  bernoulli_mock <- mockery::mock(c(FALSE, TRUE, FALSE))
  mockery::stub(schedule_outcome, 'bernoulli_multi_p', bernoulli_mock)
  schedule_outcome(
    api = api,
    target = c(5, 6, 7),
    prob_successful = c(.3, .4, .5),
    success_event = success_event,
    failure_event = failure_event
  )
  mockery::expect_args(bernoulli_mock, 1, c(.3, .4, .5))
  mockery::expect_args(
    api$schedule,
    1,
    event = success_event,
    target = 6,
    delay = 0
  )
  mockery::expect_args(
    api$schedule,
    2,
    event = failure_event,
    target = c(5, 7),
    delay = 0
  )
})

test_that('allocate_treatment can allocate for one individual', {
  api <- list(get_state = mockery::mock(c(1, 5, 7))) # mock treated
  human <- mockery::mock()

  allocated <- allocate_treatment(
    api = api,
    human = human,
    need_treatment = 4,
    treated_state = mockery::mock(),
    limit = 5
  )

  expect_equal(allocated, 4)
})

test_that('allocate_treatment can allocate when limit is exceeded', {
  api <- list(get_state = mockery::mock(c(1, 5, 7))) # mock treated
  human <- mockery::mock()

  mockery::stub(allocate_treatment, 'sample.int', mockery::mock(c(1, 3)))

  allocated <- allocate_treatment(
    api = api,
    human = human,
    need_treatment = c(2, 3, 4),
    treated_state = mockery::mock(),
    limit = 5
  )

  expect_equal(allocated, c(2, 4))
})

test_that('hospitilisation_flow_process can integrate a mix of mv and ox', {
  api <- list(
    get_variable = mockery::mock(c(1, 2, 3, 4)), # mock age
    get_parameters = mockery::mock(list(
      prob_severe_death_treatment = rep(.1, 4),
      prob_severe_death_no_treatment = rep(.2, 4),
      prob_non_severe_death_treatment = rep(.3, 4),
      prob_non_severe_death_no_treatment = rep(.4, 4)
    ))
  )

  events <- create_events()

  process <- hospitilisation_flow_process(
    discrete_age = mockery::mock(),
    human = mockery::mock(),
    states = list(),
    events = events
  )

  mockery::stub(
    process,
    'bernoulli_multi_p',
    mockery::mock(c(TRUE, TRUE, FALSE, FALSE))
  )

  mockery::stub(
    process,
    'allocate_treatment',
    mockery::mock(1, 3)
  )

  schedule_mock <- mockery::mock()
  mockery::stub(
    process,
    'schedule_outcome',
    schedule_mock
  )

  process(api, hospitalised = c(1, 2, 3, 4))

  mockery::expect_args(
    schedule_mock,
    1,
    api = api,
    target = 1,
    prob_successful = .1,
    success_event = events$imv_get_die,
    failure_event = events$imv_get_live
  )

  mockery::expect_args(
    schedule_mock,
    2,
    api = api,
    target = 2,
    prob_successful = .2,
    success_event = events$imv_not_get_die,
    failure_event = events$imv_not_get_live
  )

  mockery::expect_args(
    schedule_mock,
    3,
    api = api,
    target = 3,
    prob_successful = .3,
    success_event = events$iox_get_die,
    failure_event = events$iox_get_live
  )

  mockery::expect_args(
    schedule_mock,
    4,
    api = api,
    target = 4,
    prob_successful = .4,
    success_event = events$iox_not_get_die,
    failure_event = events$iox_not_get_live
  )
})
