# test_that("vaccination_process phase 1 step 1 working", {
#
#   variables <- create_vaccine_variables(variables = variables,pop = n,max_dose = parameters$N_phase)
#   events <- list(scheduled_dose = replicate(n = parameters$N_phase,expr = individual::TargetedEvent$new(n),simplify = FALSE))
#
#   vax_proc <- vaccination_process(parameters = parameters,variables = variables,events = events,dt = 1)
#   vax_proc(timestep = 1)
#
#   sched <- events$scheduled_dose[[1]]$get_scheduled()$copy()
#
#   expect_true(
#     all(variables$discrete_age$get_values(sched) == 17)
#   )
#   expect_equal(
#     sched$size(), parameters$vaccine_set[1]
#   )
#   expect_equal(
#     events$scheduled_dose[[2]]$get_scheduled()$size(), 0
#   )
#   expect_equal(
#     events$scheduled_dose[[3]]$get_scheduled()$size(), 0
#   )
#
# })
#
#
# test_that("vaccination_process phase 1 step 1; assigns for next step if nobody eligible", {
#
#   variables <- create_vaccine_variables(variables = variables,pop = n,max_dose = parameters$N_phase)
#   events <- list(scheduled_dose = replicate(n = parameters$N_phase,expr = individual::TargetedEvent$new(n),simplify = FALSE))
#
#   schedule_dose_vaccine(timestep = 1,variables = variables,target = which(ages==17),dose = 1)
#   update_vaccine_variables(variables = variables)
#
#   # should schedule doses for 1st dose in age group 16 (no priority doses because it's too soon)
#   vax_proc <- vaccination_process(parameters = parameters,variables = variables,events = events,dt = 1)
#   vax_proc(timestep = 2)
#
#   expect_true(
#     all(variables$discrete_age$get_values(events$scheduled_dose[[1]]$get_scheduled()) == 16)
#   )
#   expect_equal(
#     events$scheduled_dose[[1]]$get_scheduled()$size(), parameters$vaccine_set[2]
#   )
#   expect_equal(
#     events$scheduled_dose[[2]]$get_scheduled()$size(), 0
#   )
#   expect_equal(
#     events$scheduled_dose[[3]]$get_scheduled()$size(), 0
#   )
#
# })
#
#
# test_that("vaccination_process phase 1 step 1; assigns priority doses to next step", {
#
#   variables <- create_vaccine_variables(variables = variables,pop = n,max_dose = parameters$N_phase)
#   events <- list(scheduled_dose = replicate(n = parameters$N_phase,expr = individual::TargetedEvent$new(n),simplify = FALSE))
#
#   schedule_dose_vaccine(timestep = 1,variables = variables,target = which(ages==17),dose = 1)
#   update_vaccine_variables(variables = variables)
#
#   vax_proc <- vaccination_process(parameters = parameters,variables = variables,events = events,dt = 1)
#   vax_proc(timestep = parameters$dose_period[2] + 1)
#
#   expect_true(
#     all(variables$discrete_age$get_values(events$scheduled_dose[[2]]$get_scheduled()) == 17)
#   )
#   expect_equal(
#     events$scheduled_dose[[1]]$get_scheduled()$size(), 0
#   )
#   expect_equal(
#     events$scheduled_dose[[2]]$get_scheduled()$size(),  parameters$vaccine_set[parameters$dose_period[2] + 1]
#   )
#   expect_equal(
#     events$scheduled_dose[[3]]$get_scheduled()$size(), 0
#   )
#
# })
#
#
# test_that("vaccination_process phase 1 step 17; move to next phase (2) and vaccinate persons there (step 4)", {
#
#   variables <- create_vaccine_variables(variables = variables,pop = n,max_dose = parameters$N_phase)
#   events <- list(scheduled_dose = replicate(n = parameters$N_phase,expr = individual::TargetedEvent$new(n),simplify = FALSE))
#
#   # give everyone 1st dose
#   schedule_dose_vaccine(timestep = 1,variables = variables,target = full_bset,dose = 1)
#   update_vaccine_variables(variables = variables)
#   # give prioritized groups 2nd dose
#   pri_bset <- variables$discrete_age$get_index_of(set = which(parameters$next_dose_priority[1, ] > 0))
#   schedule_dose_vaccine(timestep = 2,variables = variables,target = pri_bset,dose = 2)
#   update_vaccine_variables(variables = variables)
#
#   vax_proc <- vaccination_process(parameters = parameters,variables = variables,events = events,dt = 1)
#   vax_proc(timestep = parameters$dose_period[2] + 1)
#
#   expect_equal(
#     variables$phase$value, 2
#   )
#
#   expect_true(
#     all(variables$discrete_age$get_values(events$scheduled_dose[[2]]$get_scheduled()) == 14)
#   )
#   expect_equal(
#     events$scheduled_dose[[1]]$get_scheduled()$size(), 0
#   )
#   expect_equal(
#     events$scheduled_dose[[2]]$get_scheduled()$size(),  parameters$vaccine_set[parameters$dose_period[2] + 1]
#   )
#   expect_equal(
#     events$scheduled_dose[[3]]$get_scheduled()$size(), 0
#   )
#
# })
#
#
# test_that("vaccination_process does not do anything if done with all phases", {
#
#   variables <- create_vaccine_variables(variables = variables,pop = n,max_dose = parameters$N_phase)
#   events <- list(scheduled_dose = replicate(n = parameters$N_phase,expr = individual::TargetedEvent$new(n),simplify = FALSE))
#
#   # give everyone all doses
#   schedule_dose_vaccine(timestep = 1,variables = variables,target = full_bset,dose = 1)
#   schedule_dose_vaccine(timestep = 2,variables = variables,target = full_bset,dose = 2)
#   schedule_dose_vaccine(timestep = 3,variables = variables,target = full_bset,dose = 3)
#   update_vaccine_variables(variables = variables)
#
#   vax_proc <- vaccination_process(parameters = parameters,variables = variables,events = events,dt = 1)
#   vax_proc(timestep = 30)
#
#   expect_equal(
#     variables$phase$value, 4
#   )
#   expect_equal(
#     events$scheduled_dose[[1]]$get_scheduled()$size(), 0
#   )
#   expect_equal(
#     events$scheduled_dose[[2]]$get_scheduled()$size(),0
#   )
#   expect_equal(
#     events$scheduled_dose[[3]]$get_scheduled()$size(), 0
#   )
#
# })
#
#
# test_that("vaccination_process on phase 1 step 1 advances correctly if we have an excess of doses", {
#
#   variables <- create_vaccine_variables(variables = variables,pop = n,max_dose = parameters$N_phase)
#   events <- list(scheduled_dose = replicate(n = parameters$N_phase,expr = individual::TargetedEvent$new(n),simplify = FALSE))
#
#   parameters$vaccine_set[10] <- 80*3
#   vax_proc <- vaccination_process(parameters = parameters,variables = variables,events = events,dt = 1)
#   vax_proc(timestep = 10)
#
#   vaxx_by_age <- tab_bins(variables$discrete_age$get_values(events$scheduled_dose[[1]]$get_scheduled()), 17)
#
#   expect_true(all(vaxx_by_age < 83)) # check it's giving a reasonable number to each age group
#   expect_equal(
#     events$scheduled_dose[[1]]$get_scheduled()$size(), parameters$vaccine_set[10]
#   )
#   expect_equal(
#     events$scheduled_dose[[2]]$get_scheduled()$size(), 0
#   )
#   expect_equal(
#     events$scheduled_dose[[3]]$get_scheduled()$size(), 0
#   )
#
# })
