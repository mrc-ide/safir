# test_that('prob_outcome draws correct probabilities', {
#   expect_equal(
#     prob_outcome(
#       age = c(1, 1, 1, 2, 3, 4, 5),
#       probs = c(.1, .2, .3, .4, .5, .6, .7)
#     ),
#     c(.1, .2, .3, .4, .5, .6, .7)[c(1, 1, 1, 2, 3, 4, 5)]
#   )
# })
#
# test_that('schedule_outcome can schedule both successes and failures', {
#
#   success_event <- mock_event()
#   failure_event <- mock_event()
#   bernoulli_mock <- mockery::mock(c(FALSE, TRUE, FALSE))
#   mockery::stub(schedule_outcome, 'bernoulli_multi_p', bernoulli_mock)
#   target <- individual::Bitset$new(3)
#   target <- target$insert(1:3)
#
#   schedule_outcome(
#     target,
#     prob_successful = c(.3, .4, .5),
#     success_event = success_event,
#     failure_event = failure_event
#   )
#
#   mockery::expect_args(bernoulli_mock, 1, c(.3, .4, .5))
#
#   to_succeed <- individual::filter_bitset(target, 2)
#   mockery::expect_args(
#     success_event$schedule,
#     1,
#     target = to_succeed,
#     delay = 0
#   )
#   to_fail <- individual::filter_bitset(target, c(1, 3))
#   mockery::expect_args(
#     failure_event$schedule,
#     1,
#     target = to_fail,
#     delay = 0
#   )
# })
#
# test_that('allocate_treatment can allocate when limit is not exceeded', {
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
#   variables$states <- mock_category(c("S", 'IOxGetDie', 'IOxGetLive', 'IRec'),
#                                     c("S", "S", "S", "S", "S"))
#
#   to_treat <- individual::Bitset$new(5)
#   to_treat <- to_treat$insert(1:5)
#   allocated <- allocate_treatment(
#     variables,
#     to_treat,
#     c('IOxGetDie', 'IOxGetLive', 'IRec'),
#     100
#   )
#
#   expect_equal(allocated$size(), 5)
# })
#
# test_that('allocate_treatment can allocate when limit is exceeded', {
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
#   variables$states <- mock_category(c("S", 'IOxGetDie', 'IOxGetLive', 'IRec'),
#                                     c("S", "S", "IRec", "IRec", "IRec"))
#
#   to_treat <- individual::Bitset$new(5)
#   to_treat <- to_treat$insert(1:2)
#   allocated <- allocate_treatment(
#     variables,
#     to_treat,
#     c('IOxGetDie', 'IOxGetLive', 'IRec'),
#     3
#   )
#
#   expect_equal(allocated$size(), 0)
# })
#
# test_that('hospitilisation_flow_process can integrate a mix of mv and ox', {
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
#   events <- create_events(parameters)
#   variables <- create_variables(pop, parameters)
#   attach_event_listeners(
#     variables,
#     events,
#     parameters
#   )
#
#   process <- hospitilisation_flow_process(
#     parameters,
#     variables,
#     events
#   )
#
#   variables$discrete_age <- mock_integer(rep(1, 5))
#
#   # First 2 need MV second 2 Ox
#   mockery::stub(
#     process,
#     'bernoulli_multi_p',
#     mockery::mock(c(TRUE, TRUE, FALSE, FALSE))
#   )
#
#   get_treat1 <- individual::Bitset$new(4)
#   get_treat1 <- get_treat1$insert(1)
#   get_treat2 <- individual::Bitset$new(4)
#   get_treat2 <- get_treat2$insert(4)
#   mockery::stub(
#     process,
#     'allocate_treatment',
#     mockery::mock(get_treat1, get_treat2)
#   )
#
#   schedule_mock <- mockery::mock()
#   mockery::stub(
#     process,
#     'schedule_outcome',
#     schedule_mock
#   )
#
#   hospitalised <- individual::Bitset$new(4)
#   hospitalised <- hospitalised$insert(1:4)
#   process(0, hospitalised)
#
#
#   t1 <- individual::Bitset$new(4)
#   t1 <- t1$insert(1)
#   mockery::expect_args(
#     schedule_mock,
#     1,
#     target = t1,
#     prob_successful = parameters$prob_severe_death_treatment[1],
#     success_event = events$imv_get_die,
#     failure_event = events$imv_get_live
#   )
#
#   t2 <- individual::Bitset$new(4)
#   t2 <- t2$insert(2)
#   mockery::expect_args(
#     schedule_mock,
#     2,
#     target = t2,
#     prob_successful = parameters$prob_severe_death_no_treatment[1],
#     success_event = events$imv_not_get_die,
#     failure_event = events$imv_not_get_live
#   )
#
#   t3 <- individual::Bitset$new(4)
#   t3 <- t3$insert(3)
#   mockery::expect_args(
#     schedule_mock,
#     3,
#     target = t3,
#     prob_successful = parameters$prob_non_severe_death_treatment[1],
#     success_event = events$iox_get_die,
#     failure_event = events$iox_get_live
#   )
#
#   t4 <- individual::Bitset$new(4)
#   t4 <- t4$insert(4)
#   mockery::expect_args(
#     schedule_mock,
#     4,
#     target = t4,
#     prob_successful = parameters$prob_non_severe_death_no_treatment[1],
#     success_event = events$iox_not_get_die,
#     failure_event = events$iox_not_get_live
#   )
# })
