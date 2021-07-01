test_that("vaccination_process phase 1 step 1 working", {

  variables <- create_vaccine_variables(variables = variables,pop = n,max_dose = parameters$N_phase)
  events <- list(scheduled_dose = replicate(n = parameters$N_phase,expr = individual::TargetedEvent$new(n),simplify = FALSE))

  vax_proc <- vaccination_process(parameters = parameters,variables = variables,events = events,dt = 1)
  vax_proc(timestep = 1)

  expect_true(
    all(variables$discrete_age$get_values(events$scheduled_dose[[1]]$get_scheduled()) == 17)
  )
  expect_equal(
    events$scheduled_dose[[1]]$get_scheduled()$size(), parameters$vaccine_set[1]
  )
  expect_equal(
    events$scheduled_dose[[2]]$get_scheduled()$size(), 0
  )
  expect_equal(
    events$scheduled_dose[[3]]$get_scheduled()$size(), 0
  )

})


test_that("vaccination_process phase 1 step 1; assigns for next step if nobody eligible", {

  variables <- create_vaccine_variables(variables = variables,pop = n,max_dose = parameters$N_phase)
  events <- list(scheduled_dose = replicate(n = parameters$N_phase,expr = individual::TargetedEvent$new(n),simplify = FALSE))

  variables$dose_time[[1]]$queue_update(values = 1,index = which(ages==17))
  variables$dose_time[[1]]$.update()

  # should schedule doses for 1st dose in age group 16 (no priority doses because it's too soon)
  vax_proc <- vaccination_process(parameters = parameters,variables = variables,events = events,dt = 1)
  vax_proc(timestep = 2)

  expect_true(
    all(variables$discrete_age$get_values(events$scheduled_dose[[1]]$get_scheduled()) == 16)
  )
  expect_equal(
    events$scheduled_dose[[1]]$get_scheduled()$size(), parameters$vaccine_set[2]
  )
  expect_equal(
    events$scheduled_dose[[2]]$get_scheduled()$size(), 0
  )
  expect_equal(
    events$scheduled_dose[[3]]$get_scheduled()$size(), 0
  )

})


test_that("vaccination_process phase 1 step 1; assigns priority doses to next step", {

  variables <- create_vaccine_variables(variables = variables,pop = n,max_dose = parameters$N_phase)
  events <- list(scheduled_dose = replicate(n = parameters$N_phase,expr = individual::TargetedEvent$new(n),simplify = FALSE))

  variables$dose_time[[1]]$queue_update(values = 1,index = which(ages==17))
  variables$dose_time[[1]]$.update()

  vax_proc <- vaccination_process(parameters = parameters,variables = variables,events = events,dt = 1)
  vax_proc(timestep = parameters$dose_period[2] + 1)

  expect_true(
    all(variables$discrete_age$get_values(events$scheduled_dose[[2]]$get_scheduled()) == 17)
  )
  expect_equal(
    events$scheduled_dose[[1]]$get_scheduled()$size(), 0
  )
  expect_equal(
    events$scheduled_dose[[2]]$get_scheduled()$size(),  parameters$vaccine_set[parameters$dose_period[2] + 1]
  )
  expect_equal(
    events$scheduled_dose[[3]]$get_scheduled()$size(), 0
  )

})


test_that("vaccination_process phase 1 step 17; move to next phase (2) and vaccinate persons there (step 4)", {

  variables <- create_vaccine_variables(variables = variables,pop = n,max_dose = parameters$N_phase)
  events <- list(scheduled_dose = replicate(n = parameters$N_phase,expr = individual::TargetedEvent$new(n),simplify = FALSE))

  # give everyone 1st dose
  variables$dose_time[[1]]$queue_update(values = 1)
  variables$dose_time[[1]]$.update()
  # give prioritized groups 2nd dose
  pri_bset <- variables$discrete_age$get_index_of(set = which(parameters$next_dose_priority[1, ] > 0))
  variables$dose_time[[2]]$queue_update(values = 2,index = pri_bset)
  variables$dose_time[[2]]$.update()

  vax_proc <- vaccination_process(parameters = parameters,variables = variables,events = events,dt = 1)
  vax_proc(timestep = parameters$dose_period[2] + 1)

  expect_equal(
    variables$phase$value, 2
  )

  expect_true(
    all(variables$discrete_age$get_values(events$scheduled_dose[[2]]$get_scheduled()) == 14)
  )
  expect_equal(
    events$scheduled_dose[[1]]$get_scheduled()$size(), 0
  )
  expect_equal(
    events$scheduled_dose[[2]]$get_scheduled()$size(),  parameters$vaccine_set[parameters$dose_period[2] + 1]
  )
  expect_equal(
    events$scheduled_dose[[3]]$get_scheduled()$size(), 0
  )

})
