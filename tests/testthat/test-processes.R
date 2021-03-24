test_that("test create_setup_process", {
   # create our dummy variables for create_setup_process
   # these are mock arguments to be passed to the function
   # Because create_steup_process is istelf a wrapper to a function that takes
   # an argument api, we will mock the arguments passed to the api

   states <- list(
      E = mockery::mock()
   )

   human <- mockery::mock()

   ret <- create_setup_process(
      human,
      states,
      events,
      variables
   )

   problist <- seq(1,17,1)

   api <- list(get_variable = mockery::mock(),
               get_parameters = mockery::mock(),
               get_state = mockery::mock(c(6, 7, 3, 2, 5, 8)),
               schedule = mockery::mock())

   variables <- list(discrete_age = list(mockery::mock(problist)))

   events <- list(mild_infection = mockery::mock(),
                  asymp_infection = mockery::mock(),
                  severe_infection = mockery::mock())

   parameters = mockery::mock(list(dur_E = mockery::mock(),
                                   prob_hosp = list(mockery::mock(problist)),
                                   prob_asymp = list(mockery::mock(problist))))

   mockery::stub(ret, 'r_erlang', mockery::mock(c(0.5, 0.5), c(0.4, 0.4),
                                                c(0.7, 0.7)))

   mockery::stub(ret, 'bernoulli_multi_p',
                 mockery::mock(c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE),
                               c(TRUE, TRUE, FALSE, FALSE)))

   # now that we have our api let's call our ret function
   ret(api)

   # The ret should have called api$schedule three times.

   # The first call is to schedule the severe infection for individuals 7,8
   # who schedule in 0 days
   mockery::expect_args(
      api$schedule,
      1,
      event = events$severe_infection,
      target = c(6, 7),
      delay = c(0.5, 0.5)
   )

   # The second call is to schedule the severe infection for individuals 5,6
   # who schedule in 0 days
   mockery::expect_args(
      api$schedule,
      2,
      event = events$mild_infection,
      target = c(3, 2),
      delay = c(0.4, 0.4)
   )

   # The third  call is to schedule the asymptomatic infection for individuals
   # 8 9 who schedule in 0 days
   mockery::expect_args(
      api$schedule,
      3,
      event = events$asymp_infection,
      target = c(5, 8),
      delay = c(0.7, 0.7)
   )

})

test_that("test that create_pocesses works", {

   pop <- get_population("AFG")
   pop$n <- as.integer(pop$n / 10000)

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

   variables <- create_variables(pop, psq)
   states <- create_states(psq)
   events <- create_events()
   individuals <- create_human(states,
                               variables,
                               events)

   # Check create_processes has worked correctly
   for (process in create_processes(individuals,
                                    states,
                                    events,
                                    variables,
                                    psq)) {

      expect(is.function(process) || inherits(process, 'externalptr'),
             'Process is not a function')
   }

})
