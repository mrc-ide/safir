test_that("hospitalization event scheduling works, vector probabilities", {

  dt <- 1
  N_age <- 2
  N <- 2e4
  ages <- rep(1:N_age, each = N/N_age)

  parameters <- list(
    dt = dt,
    hosp_beds = N*2,
    ICU_beds = N*2,
    population = N
  )

  ICU_states <- c('IMVGetDie', 'IMVGetLive')
  hosp_states <- c('IOxGetDie', 'IOxGetLive', 'IRec')

  variables <- list(
    discrete_age = IntegerVariable$new(ages),
    states = CategoricalVariable$new(categories = c("S", ICU_states, hosp_states), initial_values = rep("S", N))
  )

  events <- create_events(parameters = parameters)

  parameters$prob_severe <- c(0.5, 0.1)
  parameters$prob_severe_death_treatment <- c(0.9, 0.99)
  parameters$prob_severe_death_no_treatment <- c(0.9, 0.99)
  parameters$prob_non_severe_death_treatment <- c(0.9, 0.99)
  parameters$prob_non_severe_death_no_treatment <- c(0.9, 0.99)

  hosp_scheduler <- create_hospital_scheduler_listener(parameters = parameters, variables = variables, events = events)

  # all empty
  expect_true(all(vapply(X = events, FUN = function(e) {
    e$get_scheduled()$size() < 1
  }, FUN.VALUE = logical(1), USE.NAMES = FALSE)))

  # schedule the pop
  hosp_scheduler(timestep = 1, hospitalised = Bitset$new(size = N)$insert(1:N))

  # check results
  sizes <- vapply(X = events, FUN = function(e) {
    e$get_scheduled()$size()
  }, FUN.VALUE = numeric(1), USE.NAMES = TRUE)

  expect_equal(sum(sizes), N)
  expect_true(abs(sum(sizes[c("imv_get_live", "imv_get_die")]) - sum(N/N_age * parameters$prob_severe)) / sum(N/N_age * parameters$prob_severe) < 0.05)
  expect_true(abs(sum(sizes[c("iox_get_live", "iox_get_die")]) - sum(N/N_age * (1 - parameters$prob_severe))) / sum(N/N_age * parameters$prob_severe) < 0.05)
  expect_true(sizes[["imv_get_die"]]/sizes[["imv_get_live"]] > 1)
  expect_true(sizes[["iox_get_die"]]/sizes[["iox_get_live"]] > 1)
})


test_that("hospitalization event scheduling works, matrix probabilities", {

  dt <- 1
  N_age <- 2
  N <- 2e4
  ages <- rep(1:N_age, each = N/N_age)

  parameters <- list(
    dt = dt,
    hosp_beds = N*2,
    ICU_beds = N*2,
    population = N
  )

  ICU_states <- c('IMVGetDie', 'IMVGetLive')
  hosp_states <- c('IOxGetDie', 'IOxGetLive', 'IRec')

  variables <- list(
    discrete_age = IntegerVariable$new(ages),
    states = CategoricalVariable$new(categories = c("S", ICU_states, hosp_states), initial_values = rep("S", N))
  )

  events <- create_events(parameters = parameters)

  parameters$prob_severe <- matrix(c(0.5, 0.1))
  parameters$prob_severe_death_treatment <- matrix(c(0.9, 0.99))
  parameters$prob_severe_death_no_treatment <- matrix(c(0.9, 0.99))
  parameters$prob_non_severe_death_treatment <- matrix(c(0.9, 0.99))
  parameters$prob_non_severe_death_no_treatment <- matrix(c(0.9, 0.99))

  hosp_scheduler <- create_hospital_scheduler_listener(parameters = parameters, variables = variables, events = events)

  # all empty
  expect_true(all(vapply(X = events, FUN = function(e) {
    e$get_scheduled()$size() < 1
  }, FUN.VALUE = logical(1), USE.NAMES = FALSE)))

  # schedule the pop
  hosp_scheduler(timestep = 1, hospitalised = Bitset$new(size = N)$insert(1:N))

  # check results
  sizes <- vapply(X = events, FUN = function(e) {
    e$get_scheduled()$size()
  }, FUN.VALUE = numeric(1), USE.NAMES = TRUE)

  expect_equal(sum(sizes), N)
  expect_true(abs(sum(sizes[c("imv_get_live", "imv_get_die")]) - sum(N/N_age * parameters$prob_severe)) / sum(N/N_age * parameters$prob_severe) < 0.05)
  expect_true(abs(sum(sizes[c("iox_get_live", "iox_get_die")]) - sum(N/N_age * (1 - parameters$prob_severe))) / sum(N/N_age * parameters$prob_severe) < 0.05)
  expect_true(sizes[["imv_get_die"]]/sizes[["imv_get_live"]] > 1)
  expect_true(sizes[["iox_get_die"]]/sizes[["iox_get_live"]] > 1)
})

