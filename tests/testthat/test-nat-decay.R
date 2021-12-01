test_that("test NAT decay from vaccination works on dt = 1", {

  library(individual)

  n <- 10
  tmax <- 800
  parameters <- get_vaccine_ab_titre_parameters(vaccine = "Pfizer")

  compare_nt <- draw_nt(parameters = parameters, n = n, tmax = tmax)

  nt <- compare_nt$nt
  z1 <- compare_nt$z1

  # safir
  parameters$population <- n
  parameters$N_phase <- 1

  variables <- create_vaccine_variables(variables = list(),parameters = parameters)

  dt <- 1
  timesteps <- tmax/dt

  schedule_dose_vaccine(timestep = 0,variables = variables,target = Bitset$new(n)$insert(1:n),dose = 1,parameters = parameters)
  variables$ab_titre$queue_update(values = log(10^z1)) # make sure using the same set of RVs
  variables$ab_titre$.update()
  variables$dose_num$.update()
  variables$dose_time$.update()

  ab_titre <- vaccine_ab_titre_process(parameters = parameters,variables = variables, dt = dt)

  safir_out <- matrix(data = NaN,nrow = timesteps + 1,ncol = n)
  safir_out[1, ] <- variables$ab_titre$get_values()

  for (t in 1:timesteps) {
    ab_titre(timestep = t)
    variables$ab_titre$.update()
    safir_out[t + 1L, ] <- variables$ab_titre$get_values()
  }

  # test Ab time series the same
  expect_equal(safir_out[nrow(safir_out), ], nt[nrow(nt), ])

})


test_that("test NAT decay from vaccination works on dt = 0.5", {

  library(individual)

  n <- 10
  tmax <- 800
  parameters <- get_vaccine_ab_titre_parameters(vaccine = "Pfizer")

  compare_nt <- draw_nt(parameters = parameters, n = n, tmax = tmax)

  nt <- compare_nt$nt
  z1 <- compare_nt$z1

  # safir
  parameters$population <- n
  parameters$N_phase <- 1

  variables <- create_vaccine_variables(variables = list(),parameters = parameters)

  dt <- 0.5
  timesteps <- tmax/dt

  schedule_dose_vaccine(timestep = 0,variables = variables,target = Bitset$new(n)$insert(1:n),dose = 1,parameters = parameters)
  variables$ab_titre$queue_update(values = log(10^z1)) # make sure using the same set of RVs
  variables$ab_titre$.update()
  variables$dose_num$.update()
  variables$dose_time$.update()

  ab_titre <- vaccine_ab_titre_process(parameters = parameters,variables = variables, dt = dt)

  safir_out <- matrix(data = NaN,nrow = timesteps + 1,ncol = n)
  safir_out[1, ] <- variables$ab_titre$get_values()

  for (t in 1:timesteps) {
    ab_titre(timestep = t)
    variables$ab_titre$.update()
    safir_out[t + 1L, ] <- variables$ab_titre$get_values()
  }

  # test Ab time series the same
  expect_equal(safir_out[nrow(safir_out), ], nt[nrow(nt), ])

})


test_that("test NAT decay from infection works on dt = 0.5", {

  library(individual)

  n <- 10
  tmax <- 800
  parameters <- get_vaccine_ab_titre_parameters(vaccine = "Pfizer")

  compare_nt <- draw_nt(parameters = parameters, n = n, tmax = tmax)

  nt <- compare_nt$nt
  z1 <- compare_nt$z1

  # safir
  parameters$population <- n
  parameters$N_phase <- 1

  variables <- create_vaccine_variables(variables = list(),parameters = parameters)
  variables <- create_natural_immunity_variables(variables = variables, parameters = parameters)

  dt <- 0.5
  timesteps <- tmax/dt

  target <- Bitset$new(n)$insert(1:n)
  variables$inf_num$queue_update(values = 1, index = target)
  variables$inf_time$queue_update(values = 0, index = target)
  variables$ab_titre$queue_update(values = log(10^z1)) # make sure using the same set of RVs
  variables$ab_titre$.update()
  variables$inf_num$.update()
  variables$inf_time$.update()

  ab_titre <- natural_immunity_ab_titre_process(parameters = parameters,variables = variables, dt = dt)

  safir_out <- matrix(data = NaN,nrow = timesteps + 1,ncol = n)
  safir_out[1, ] <- variables$ab_titre$get_values()

  for (t in 1:timesteps) {
    ab_titre(timestep = t)
    variables$ab_titre$.update()
    safir_out[t + 1L, ] <-   variables$ab_titre$get_values()
  }

  # test Ab time series the same
  expect_equal(safir_out[nrow(safir_out), ], nt[nrow(nt), ])

})
