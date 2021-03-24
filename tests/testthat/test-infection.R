
test_that("infection_process does not schedule infections when no-one is infected", {
  human <- mockery::mock()

  exposure = mockery::mock()
  discrete_age = mockery::mock()
  age = mockery::mock()

  states <- list(
    IMild = mockery::mock(),
    ICase = mockery::mock(),
    S = mockery::mock()
  )

  api <- list(
    get_state = mockery::mock(NULL),
    schedule = mockery::mock(),
    get_parameters = mockery::mock()
  )

  process <- infection_process(human, states, discrete_age, exposure)

  process(api)

  mockery::expect_args(
    api$get_state,
    1,
    human,
    states$IMild,
    states$IAsymp,
    states$ICase
  )

  mockery::expect_called(api$schedule, 0)
})

test_that("infection_process does not schedule when no-one infects by chance", {
  human <- mockery::mock()
  states <- list(
    IMild = mockery::mock(),
    ICase = mockery::mock(),
    S = mockery::mock()
  )
  exposure = mockery::mock()
  discrete_age = mockery::mock()
  age = mockery::mock()
  contact_matrix <- mockery::mock()

  process <- infection_process(human, states, discrete_age, exposure,
                               contact_matrix)

  api <- list(
    get_state = mockery::mock(
      c(1, 2, 3, 4), #infected
      c(5, 6, 7) #susceptible
    ),
    get_variable = mockery::mock(c(1, 2, 3, 3), c(2, 3, 1)), # age Inf, Susc
    get_parameters = mockery::mock(list(
      N_age = 3,
      beta = .3
    )),
    get_timestep = mockery::mock(1),
    schedule = mockery::mock()
  )

  mockery::stub(process, 'bernoulli_multi_p', mockery::mock(rep(FALSE, 3)))
  mockery::stub(
    process,
    'get_contact_matrix',
    mockery::mock(array(c(
      .1, .2, .3,
      .1, .2, .3,
      .1, .2, .3
    ), dim=c(3, 3)))
  )

  process(api)

  mockery::expect_called(api$schedule, 0)
})

test_that("infection_process gives correct lambda for infections", {
  human <- mockery::mock()
  states <- list(
    IMild = mockery::mock(),
    ICase = mockery::mock(),
    S = mockery::mock()
  )
  exposure = mockery::mock()
  discrete_age = mockery::mock()
  age = mockery::mock()
  contact_matrix <- mockery::mock()

  process <- infection_process(human, states, discrete_age, exposure,
                               contact_matrix)

  api <- list(
    get_state = mockery::mock(
      c(1, 2, 3, 4), #infected
      c(5, 6, 7) #susceptible
    ),
    get_variable = mockery::mock(c(1, 2, 3, 3), c(2, 3, 1)), # age Inf, Susc
    get_parameters = mockery::mock(list(
      N_age = 3,
      beta = .3
    )),
    get_timestep = mockery::mock(1),
    schedule = mockery::mock()
  )

  bernoulli_mock <- mockery::mock(c(FALSE, TRUE, TRUE))
  mockery::stub(process, 'bernoulli_multi_p', bernoulli_mock)
  mockery::stub(
    process,
    'get_contact_matrix',
    mockery::mock(array(c(
      .1, .2, .3,
      .1, .2, .3,
      .1, .2, .3
    ), dim=c(3, 3)))
  )

  process(api)

  expect_equal(
    mockery::mock_args(bernoulli_mock)[[1]][[1]],
    c(0.2133721, 0.3023237, 0.1130796),
    tolerance=1e-6
  )

  mockery::expect_args(
    api$schedule,
    1,
    event = exposure,
    target = c(6, 7),
    delay = 0
  )

})
