
pop <- get_population("ATG")

pop$n <- as.integer(pop$n / 100)

psq <- get_parameters(
  iso3c = "ATG",
  population = pop$n,
  R0 = 2,
  time_period = 100,
  tt_contact_matrix = 0
)

events <- create_events(psq)
variables <- create_variables(pop, psq)
attach_event_listeners(
  variables,
  events,
  psq
)


test_that("create_event_based_processes assigns a listener to each event", {
  for (event in events) {
    expect_gt(length(event$.listeners), 0)
  }
})

test_that("test create_infection_update_listener", {

  variables$states <- mock_category(c("S", "IMild"), rep("S", 5))

  ret <- create_infection_update_listener(
    variables,
    "IMild"
  )

  to_move <- individual::Bitset$new(5)
  to_move <- to_move$insert(1)

  ret(0, to_move)
  mockery::expect_args(variables$states$queue_update, 1, "IMild", to_move)

})

test_that("test create_progression_listener", {

  # Without shift
  events$exposure <- mock_event()
  ret <- create_progression_listener(events$exposure, duration = 1, shift = 0, func = r_erlang)
  # Stub so that delays = 1:4
  mockery::stub(r_erlang, "rgamma", 1:4)

  to_move <- individual::Bitset$new(4)
  to_move <- to_move$insert(1:4)
  ret(1, to_move)

  mockery::expect_args(events$exposure$schedule, 1, to_move, 1:4)

  # With shift
  events$exposure <- mock_event()
  ret <- create_progression_listener(events$exposure, duration = 1, shift = 1, func = r_erlang)
  mockery::stub(r_erlang, "rgamma", 1:4)

  to_move <- individual::Bitset$new(4)
  to_move <- to_move$insert(1:4)
  ret(1, to_move)

  mockery::expect_args(events$exposure$schedule, 1, to_move, 1:4 + 1)

})

test_that("test create_exposure_update_listener", {

  events$severe_infection <- mock_event()
  events$mild_infection <- mock_event()
  events$asymp_infection <- mock_event()

  ret <- create_exposure_update_listener(
    events,
    variables,
    psq
  )
  # Stub delays
  mockery::stub(ret, "r_erlang", mockery::mock(c(0.5, 0.5), c(0.4, 0.4), c(0.7, 0.7)))
  # Stub tree probs
  mockery::stub(ret, 'bernoulli_multi_p',
                mockery::mock(c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE),
                              c(TRUE, TRUE, FALSE, FALSE)))


  to_move <- individual::Bitset$new(6)
  to_move <- to_move$insert(1:6)

  ret(0, to_move)

  # The process should have called schedule three times

  # First call is to schedule the severe infection for individuals 1 & 2
  # who schedule in 0.5 + 1 days
  have_move1 <- individual::Bitset$new(6)
  have_move1 <- have_move1$insert(1:2)
  mockery::expect_args(
    events$severe_infection$schedule,
    1,
    have_move1,
    c(1.5, 1.5)
  )

  # The second call is to schedule the asymptomatic infection for individuals
  # 5 & 6 who schedule in 0.4 + 1 days
  have_move3 <- individual::Bitset$new(6)
  have_move3 <- have_move3$insert(5:6)
  mockery::expect_args(
    events$asymp_infection$schedule,
    1,
    have_move3,
    c(1.4, 1.4)
  )

  # The thirs call is to schedule the mild infection for individuals 3 & 4
  # who schedule in 0.7 + 1 days
  have_move2 <- individual::Bitset$new(6)
  have_move2 <- have_move2$insert(3:4)
  mockery::expect_args(
    events$mild_infection$schedule,
    1,
    have_move2,
    c(1.7, 1.7)
  )
})


