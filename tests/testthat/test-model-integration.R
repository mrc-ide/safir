test_that("run_simulation can parameterise and run an Afghan model for 10 days", {

  R0 <- 2
  time_period <- 10
  tt_contact_matrix <- 0

  iso3c <- "AFG"
  pop <- get_population(iso3c)
  pop$n <- as.integer(pop$n / 1000)

  psq <- get_parameters(
    population = pop$n,
    contact_matrix_set = squire::get_mixing_matrix(iso3c = iso3c),
    iso3c = iso3c,
    R0 = R0,
    time_period = time_period
  )

  output <- run_simulation(pop, psq)

  expect_equal(length(output$timestep), 10)
  expect_equal(nrow(output), 10)
  expect_equal(length(output$E_count), 10)
  expect_equal(length(output$IMild_count), 10)
  expect_equal(length(output$IMVNotGetLive_count), 10)

})

test_that("run_simulation produces a sensible dataframe", {

  # simple overrides
  overrides <- list("pop" = list(1,3),
                    "parameters" = list(1,3))

  # get our run function that we will do some mocking to
  run <- run_simulation_replicate

  # mock the function to just return a data frame made up of our overrides
  # to check they are being grabbed correctly
  mockery::stub(
    run,
    'run_simulation',
    function(pop, parameters, max_age) {
      data.frame("timestep" = 1,
                 "pop" = pop,
                 "parameters" = parameters
                 )
    }
  )

  # run our mocked function
  out <- run(repetitions = 2, overrides = overrides)

  # Now we can check the arguments are correctly filled from overrides
  expect_equal(out$repetition[1], 1)
  expect_equal(out$timestep[2], 1)
  expect_equal(out$pop[2], 3)

  expected <- data.frame(
    timestep = c(1, 1),
    pop = c(1, 3),
    parameters = c(1, 3),
    repetition= c(1, 2)
  )

  expect_equal(out, expected)

})

test_that("run_simulation is called with the correct arguments", {

  overrides <- list("pop" = list(1,3),
                    "parameters" = list(1,3))
  # get our run function that we will do some mocking to
  run <- run_simulation_replicate
  # mock the function as just a plain mock
  run_simulation_mock <- mockery::mock(cycle = TRUE)
  mockery::stub(
    run,
    'run_simulation',
    run_simulation_mock
  )

  # run our mocked function
  out <- run(repetitions = 2, overrides = overrides)

  mockery::expect_args(
    run_simulation_mock,
    1,
    pop = overrides$pop[[1]],
    parameters = overrides$parameters[[1]]
  )

  mockery::expect_args(
    run_simulation_mock,
    2,
    pop = overrides$pop[[2]],
    parameters = overrides$parameters[[2]]
  )

})

test_that("run 2 models with run_simulation sequentially on real data", {

  R0 <- 2
  time_period <- 10
  tt_contact_matrix <- 0

  iso3c <- "ATG"
  pop <- get_population(iso3c)
  pop$n <- as.integer(pop$n / 100)

  psq <- get_parameters(
    population = pop$n,
    contact_matrix_set = squire::get_mixing_matrix(iso3c = iso3c),
    iso3c = iso3c,
    R0 = R0,
    time_period = time_period
  )

  iso3c <- "AFG"
  pop2 <- get_population(iso3c)
  pop2$n <- as.integer(pop2$n / 100)

  psq2 <- get_parameters(
    population = pop2$n,
    contact_matrix_set = squire::get_mixing_matrix(iso3c = iso3c),
    iso3c = iso3c,
    R0 = R0,
    time_period = time_period
  )

  repetitions <- 2

  overrides <- list(pop = list(pop, pop2), parameters = list(psq, psq2))

  dfs <- run_simulation_replicate(
    repetitions = repetitions,
    overrides = overrides,
    parallel = FALSE
  )

  expect_equal(length(dfs), 18)
  expect_equal(length(dfs$timestep), 20)
  expect_equal(dfs$S_count[1], 952)
  expect_equal(dfs$E_count[1], 20)
  expect_equal(dfs$repetition[1], 1)
  expect_equal(dfs$timestep[1], 1)
  expect_equal(dfs$timestep[2], 2)


  # and with parallel
  options("mc.cores" = 2)
  dfs <- run_simulation_replicate(
    repetitions = repetitions,
    overrides = overrides,
    parallel = TRUE
  )

  expect_equal(length(dfs), 18)
  expect_equal(length(dfs$timestep), 20)
  expect_equal(dfs$S_count[1], 952)
  expect_equal(dfs$E_count[1], 20)
  expect_equal(dfs$repetition[1], 1)
  expect_equal(dfs$timestep[1], 1)
  expect_equal(dfs$timestep[2], 2)


})
