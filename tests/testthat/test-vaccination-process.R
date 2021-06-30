

test_that("vaccination_process working correctly", {

  parameters <- list()
  parameters$vaccine_coverage_mat <- nimue::strategy_matrix(strategy = "Elderly",max_coverage = 0.8)
  parameters$N_age <- ncol(parameters$vaccine_coverage_mat)
  parameters$N_prioritisation_steps <- nrow(parameters$vaccine_coverage_mat)
  parameters$N_phase <- 3
  parameters$dose_period <- c(NaN, 14, 7)
  parameters$vaccine_set <- rep(80,100)

  parameters$next_dose_priority <- matrix(0, nrow = parameters$N_phase - 1, ncol = parameters$N_age)
  parameters$next_dose_priority[1, 15:17] <- 1
  parameters$next_dose_priority[2, 10:14] <- 1

  n <- 17 * 100
  ages <- rep(1:17, each = 100)

  variables <- list()
  variables$discrete_age <- IntegerVariable$new(ages)

  variables <- create_vaccine_variables(variables = variables,pop = n,max_dose = parameters$N_phase)

  events <- list(scheduled_dose = replicate(n = parameters$N_phase,expr = individual::TargetedEvent$new(n),simplify = FALSE))

  # test for phase 1 step 1
  vax_proc <- vaccination_process(parameters = parameters,variables = variables,events = events,dt = 1)
  vax_proc(timestep = 1)

  events$scheduled_dose[[1]]$get_scheduled()

  expect_true(
    all(variables$discrete_age$get_values(events$scheduled_dose[[1]]$get_scheduled()) == 17)
  )
})
