test_that("test that create_human works", {

  # Get population information
  pop <- get_population("AFG")
  pop$n <- as.integer(pop$n/1000)
  R0 <- 2
  time_period <- 1000
  tt_contact_matrix <- 0

  psq <- get_parameters(
    iso3c = "AFG",
    population = pop$n,
    R0 = R0,
    time_period = time_period,
    tt_contact_matrix = tt_contact_matrix
  )

  states <- create_states(psq)
  events <- create_events()
  variables <- create_variables(pop = pop, psq)

  indiv <- create_human(states, variables, events)

  # Create test human
  human <- individual::Individual$new(
    "human",
    states = list(
      states$S,
      states$E,
      states$IMild,
      states$IAsymp,
      states$ICase,
      states$IOxGetLive,
      states$IOxGetDie,
      states$IOxNotGetLive,
      states$IOxNotGetDie,
      states$IMVGetLive,
      states$IMVGetDie,
      states$IMVNotGetLive,
      states$IMVNotGetDie,
      states$IRec,
      states$R,
      states$D),
    variables = variables,
    events = events
  )

  expect_equal(indiv$name, human$name)
  expect_equal(length(indiv$states), length(human$states))

})

test_that("test Create_states with S for 1st age group", {

  pop <- get_population("AFG")
  pop$n <- as.integer(pop$n/1000)
  R0 <- 2
  time_period <- 1000
  tt_contact_matrix <- 0

  psq <- get_parameters(
    iso3c = "AFG",
    population = pop$n,
    R0 = R0,
    time_period = time_period,
    tt_contact_matrix = tt_contact_matrix
  )

  Snew <- individual::State$new("S", sum(psq$S_0))

  states <- create_states(psq)

  expect_equal(Snew$initial_size, states[[1]]$initial_size[1])

})
