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
