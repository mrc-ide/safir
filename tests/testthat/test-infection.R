#
# test_that("infection_process does not schedule infections when no-one is infected", {
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
#   variables$states <- mock_category(c("S", "IMild", "IAsymp", "ICase"), rep("S", 5))
#   events$exposure <- mock_event()
#
#   process <- infection_process(parameters, variables, events)
#   process(0)
#
#   mockery::expect_called(events$exposure$schedule, 0)
# })
#
# test_that("infection_process does not schedule when no-one infects by chance", {
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
#   variables$states <- mock_category(c("S", "IMild", "IAsymp", "ICase"), c(rep("S", 3), "ICase", "ICase"))
#   events$exposure <- mock_event()
#
#   process <- infection_process(parameters, variables, events)
#
#   mockery::stub(process, 'bernoulli_multi_p', mockery::mock(rep(FALSE, 3)))
#
#   process(1)
#   mockery::expect_called(events$exposure$schedule, 0)
# })
#
#
# test_that("infection_process gives correct lambda for infections", {
#
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
#   parameters$N_age = 3
#   parameters$beta_set = 0.3
#
#   events <- create_events(parameters)
#   variables <- create_variables(pop, parameters)
#   attach_event_listeners(
#     variables,
#     events,
#     parameters
#   )
#   variables$states <- mock_category(c("S", "IMild", "IAsymp", "ICase"), c(rep("ICase", 4), rep("S", 3)))
#   variables$discrete_age <- mock_integer(c(1, 2, 3, 3, 2, 3, 1))
#   events$exposure <- mock_event()
#
#   process <- infection_process(parameters, variables, events)
#
#   bernoulli_mock <- mockery::mock(c(FALSE, TRUE, TRUE))
#   mockery::stub(process, 'bernoulli_multi_p', bernoulli_mock)
#   mockery::stub(
#     process,
#     'get_contact_matrix',
#     mockery::mock(array(c(
#       .1, .2, .3,
#       .1, .2, .3,
#       .1, .2, .3
#     ), dim=c(3, 3)))
#   )
#
#   process(1)
#
#   expect_equal(
#     mockery::mock_args(bernoulli_mock)[[1]][[1]],
#     c(0.2133721, 0.3023237, 0.1130796),
#     tolerance=1e-6
#   )
#
#   to_infect <- individual::Bitset$new(7)
#   to_infect <- to_infect$insert(6:7)
#   mockery::expect_args(
#     events$exposure$schedule,
#     1,
#     target = to_infect,
#     delay = 0
#   )
# })
