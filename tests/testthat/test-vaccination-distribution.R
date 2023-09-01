test_that("coverage and get_proportion_vaccinated are giving the same results", {

  dose_times <- list(matrix(c(1, 2, NA, 2, 3, NA), nrow = 3),
                     matrix(c(NA, 3, 4, NA, NA, NA), nrow = 3),
                     matrix(c(1, 2, 2, NA, NA, NA), nrow = 3))

  dose_1 <- unlist(lapply(dose_times,function(x){x[,1]}))
  dose_1[which(is.na(dose_1))] <- -1

  dose_2 <- unlist(lapply(dose_times,function(x){x[,2]}))
  dose_2[which(is.na(dose_2))] <- -1

  n <- length(dose_1)
  dose_num <- ifelse(dose_1 == -1, 0, 1)
  dose_num[which(dose_2 > -1)] <- 2

  parameters <- list(population=vapply(dose_times,nrow,integer(1)), N_phase = 2, correlated = FALSE, N_age = 3)

  variables <- list()
  variables$discrete_age <- IntegerVariable$new(rep(seq_len(length(dose_times)),times=vapply(dose_times,nrow,integer(1))))
  variables <- create_vaccine_variables(variables = variables, parameters = parameters)
  initialize_vaccine_variables(variables = variables,dose_time_init = list(dose_1,dose_2),dose_num_init = dose_num)

  events <- list(scheduled_dose = replicate(n = 2,expr = {TargetedEvent$new(sum(vapply(dose_times,nrow,integer(1))))}))

  cov_safir <- get_current_coverage(variables = variables,events = events,dose = 1,parameters = parameters)
  cov_safir <- vapply(cov_safir,function(b){b$size()},numeric(1))

  expect_equal(
    cov_safir,
    c(2,2,3)
  )

  cov_safir <- get_current_coverage(variables = variables,events = events,dose = 2,parameters = parameters)
  cov_safir <- vapply(cov_safir,function(b){b$size()},numeric(1))

  expect_equal(
    cov_safir,
    c(2,0,0)
  )
})


test_that('eligable_for_second and eligible_for_dose_vaccine give equivalent results with same input, daily time step', {

  dose_times <- list(matrix(c(1, 2, NA, NA, 3, NA), nrow = 3),
                     matrix(c(NA, 1, 1, NA, 3, NA), nrow = 3),
                     matrix(c(1, 2, 2, NA, NA, NA), nrow = 3))

  dose_1 <- unlist(lapply(dose_times,function(x){x[,1]}))
  dose_1[which(is.na(dose_1))] <- -1

  dose_2 <- unlist(lapply(dose_times,function(x){x[,2]}))
  dose_2[which(is.na(dose_2))] <- -1

  n <- length(dose_1)
  dose_num <- ifelse(dose_1 == -1, 0, 1)
  dose_num[which(dose_2 > -1)] <- 2

  ages <- rep(seq_len(length(dose_times)),times=vapply(dose_times,nrow,integer(1)))

  parameters <- list(population=vapply(dose_times,nrow,integer(1)), N_phase = 2, correlated = FALSE, N_age = 3)

  variables <- list()
  variables$discrete_age <- IntegerVariable$new(ages)
  variables <- create_vaccine_variables(variables = variables,parameters = parameters)
  initialize_vaccine_variables(variables = variables,dose_time_init = list(dose_1,dose_2),dose_num_init = dose_num)

  events <- list(scheduled_dose = replicate(n = 2,expr = {TargetedEvent$new(n)}))

  # daily time step ------------------------------------------------------------
  dt <- 1

  t <- 1
  parameters <- list(
    dose_period = c(NaN, 14),
    N_phase = 2,
    N_age = 3,
    population = tab_bins(ages,3)
  )

  cov <- get_current_coverage(variables = variables,events = events,dose = 2,parameters = parameters)
  eligible <- get_current_eligible_from_coverage(timestep = t,dt = dt,coverage = cov,variables = variables,dose = 2,parameters = parameters)

  eligible_nimue <- nimue_eligable_for_second(dose_times, t, parameters$dose_period[2])

  expect_equal(
    vapply(eligible_nimue, function(x){sum(x)}, numeric(1)),
    vapply(eligible, function(b){b$size()}, numeric(1))
  )

  t <- 1
  parameters$dose_period[2] <- 0

  cov <- get_current_coverage(variables = variables,events = events,dose = 2,parameters = parameters)
  eligible <- get_current_eligible_from_coverage(timestep = t,dt = dt,coverage = cov,variables = variables,dose = 2,parameters = parameters)

  eligible_nimue <- nimue_eligable_for_second(dose_times, t, parameters$dose_period[2])


  vapply(cov,function(b){b$size()},numeric(1))

  expect_equal(
    vapply(eligible_nimue,function(x){sum(x)},numeric(1)),
    vapply(eligible,function(b){b$size()},numeric(1))
  )

  t <- 200
  parameters$dose_period[2] <- 14

  cov <- get_current_coverage(variables = variables,events = events,dose = 2,parameters = parameters)
  eligible <- get_current_eligible_from_coverage(timestep = t,dt = dt,coverage = cov,variables = variables,dose = 2,parameters = parameters)

  eligible_nimue <- nimue_eligable_for_second(dose_times, t, parameters$dose_period[2])

  expect_equal(
    vapply(eligible_nimue,function(x){sum(x)},numeric(1)),
    vapply(eligible,function(b){b$size()},numeric(1))
  )

})


test_that('eligable_for_second and eligible_for_dose_vaccine give equivalent results with same input, sub daily time step', {

  dose_times <- list(matrix(c(1, 2, NA, NA, 3, NA), nrow = 3),
                     matrix(c(NA, 1, 1, NA, 3, NA), nrow = 3),
                     matrix(c(1, 2, 2, NA, NA, NA), nrow = 3))

  # sub-daily time step ------------------------------------------------------------
  dt <- 0.2

  dose_1 <- unlist(lapply(dose_times,function(x){x[,1]}))
  dose_1[!(dose_1 %in% c(1,-1))] <- dose_1[!(dose_1 %in% c(1,-1))] / dt
  dose_1[which(is.na(dose_1))] <- -1

  dose_2 <- unlist(lapply(dose_times,function(x){x[,2]}))
  dose_2[!(dose_2 %in% c(1,-1))] <- dose_2[!(dose_2 %in% c(1,-1))] / dt
  dose_2[which(is.na(dose_2))] <- -1

  n <- length(dose_1)
  dose_num <- ifelse(dose_1 == -1, 0, 1)
  dose_num[which(dose_2 > -1)] <- 2

  ages <- rep(seq_len(length(dose_times)),times=vapply(dose_times,nrow,integer(1)))

  parameters <- list(
    dose_period = c(NaN, 14),
    N_phase = 2,
    N_age = 3,
    population = tab_bins(ages,3),
    correlated = FALSE
  )

  variables <- list()
  variables$discrete_age <- IntegerVariable$new(ages)
  variables <- create_vaccine_variables(variables = variables,parameters = parameters)
  initialize_vaccine_variables(variables = variables,dose_time_init = list(dose_1,dose_2),dose_num_init = dose_num)

  events <- list(scheduled_dose = replicate(n = 2,expr = {TargetedEvent$new(n)}))

  t <- 1 # current day
  tdt <- t/dt # current time step

  cov <- get_current_coverage(variables = variables,events = events,dose = 2,parameters = parameters)
  eligible <- get_current_eligible_from_coverage(timestep = t,dt = dt,coverage = cov,variables = variables,dose = 2,parameters = parameters)

  eligible_nimue <- nimue_eligable_for_second(dose_times, t, parameters$dose_period[2])

  expect_equal(
    vapply(eligible_nimue,function(x){sum(x)},numeric(1)),
    vapply(eligible,function(b){b$size()},numeric(1))
  )

  t <- 1
  tdt <- t/dt
  parameters$dose_period[2] <- 0

  cov <- get_current_coverage(variables = variables,events = events,dose = 2,parameters = parameters)
  eligible <- get_current_eligible_from_coverage(timestep = t,dt = dt,coverage = cov,variables = variables,dose = 2,parameters = parameters)

  eligible_nimue <- nimue_eligable_for_second(dose_times, t, parameters$dose_period[2])

  expect_equal(
    vapply(eligible_nimue,function(x){sum(x)},numeric(1)),
    vapply(eligible,function(b){b$size()},numeric(1))
  )

  t <- 200
  tdt <- t/dt
  parameters$dose_period[2] <- 14

  cov <- get_current_coverage(variables = variables,events = events,dose = 2,parameters = parameters)
  eligible <- get_current_eligible_from_coverage(timestep = t,dt = dt,coverage = cov,variables = variables,dose = 2,parameters = parameters)

  eligible_nimue <- nimue_eligable_for_second(dose_times, t, parameters$dose_period[2])

  expect_equal(
    vapply(eligible_nimue,function(x){sum(x)},numeric(1)),
    vapply(eligible,function(b){b$size()},numeric(1))
  )

})


test_that('eligable_for_second and age_group_eligible_for_dose_vaccine give equivalent results with same input', {

  dose_times <- list(matrix(c(1, 2, NA, 4, 3, NA), nrow = 3),
                     matrix(c(NA, 3, 4, NA, NA, NA), nrow = 3),
                     matrix(c(1, 2, 2, NA, NA, NA), nrow = 3))

  dose_1 <- unlist(lapply(dose_times,function(x){x[,1]}))
  dose_1[which(is.na(dose_1))] <- -1

  dose_2 <- unlist(lapply(dose_times,function(x){x[,2]}))
  dose_2[which(is.na(dose_2))] <- -1

  n <- length(dose_1)
  dose_num <- ifelse(dose_1 == -1, 0, 1)
  dose_num[which(dose_2 > -1)] <- 2

  ages <- rep(seq_len(length(dose_times)),times=vapply(dose_times,nrow,integer(1)))

  parameters <- list(
    dose_period = c(NaN, 14),
    N_age = 3,
    N_phase = 2,
    correlated = FALSE,
    population = tab_bins(ages,3)
  )

  variables <- list()
  variables$discrete_age <- IntegerVariable$new(ages)
  variables <- create_vaccine_variables(variables = variables,parameters = parameters)
  initialize_vaccine_variables(variables = variables,dose_time_init = list(dose_1,dose_2),dose_num_init = dose_num)

  events <- list(scheduled_dose = replicate(n = 2,expr = {TargetedEvent$new(n)}))

  t <- 1
  dt <- 1

  cov <- get_current_coverage(variables = variables,events = events,dose = 2,parameters = parameters)
  eligible <- get_current_eligible_from_coverage(timestep = t,dt = dt,coverage = cov,variables = variables,dose = 2,parameters = parameters)

  eligible_nimue <- nimue_eligable_for_second(dose_times, t, parameters$dose_period[2])

  expect_equal(
    vapply(eligible_nimue,function(x){sum(x)},numeric(1)),
    vapply(eligible,function(b){b$size()},numeric(1))
  )

  t <- 1
  parameters$dose_period[2] <- 0

  cov <- get_current_coverage(variables = variables,events = events,dose = 2,parameters = parameters)
  eligible <- get_current_eligible_from_coverage(timestep = t,dt = dt,coverage = cov,variables = variables,dose = 2,parameters = parameters)

  eligible_nimue <- nimue_eligable_for_second(dose_times, t, parameters$dose_period[2])

  expect_equal(
    vapply(eligible_nimue,function(x){sum(x)},numeric(1)),
    vapply(eligible,function(b){b$size()},numeric(1))
  )

  t <- 200
  parameters$dose_period[2] <- 14

  cov <- get_current_coverage(variables = variables,events = events,dose = 2,parameters = parameters)
  eligible <- get_current_eligible_from_coverage(timestep = t,dt = dt,coverage = cov,variables = variables,dose = 2,parameters = parameters)

  eligible_nimue <- nimue_eligable_for_second(dose_times, t, parameters$dose_period[2])

  expect_equal(
    vapply(eligible_nimue,function(x){sum(x)},numeric(1)),
    vapply(eligible,function(b){b$size()},numeric(1))
  )

})


test_that("target_pop is giving the same results as nimue", {

  dose_times <- list(matrix(c(1, 2, NA, NA, 3, NA), nrow = 3),
                     matrix(c(NA, 3, 4, NA, NA, NA), nrow = 3),
                     matrix(c(1, 2, 2, NA, NA, NA), nrow = 3))

  dose_1 <- unlist(lapply(dose_times,function(x){x[,1]}))
  dose_1[which(is.na(dose_1))] <- -1

  dose_2 <- unlist(lapply(dose_times,function(x){x[,2]}))
  dose_2[which(is.na(dose_2))] <- -1

  n <- length(dose_1)
  dose_num <- ifelse(dose_1 == -1, 0, 1)
  dose_num[which(dose_2 > -1)] <- 2

  ages <- rep(seq_len(length(dose_times)),times=vapply(dose_times,nrow,integer(1)))

  parameters <- list(population=vapply(dose_times,nrow,integer(1)), N_phase = 2, correlated = FALSE)

  variables <- list()
  variables$discrete_age <- IntegerVariable$new(ages)
  variables <- create_vaccine_variables(variables = variables,parameters = parameters)
  initialize_vaccine_variables(variables = variables,dose_time_init = list(dose_1,dose_2),dose_num_init = dose_num)

  events <- list(scheduled_dose = replicate(n = 2,expr = {TargetedEvent$new(n)}))

  parameters <- list(
    N_age = 3,
    dose_period = c(NaN, 14, NaN),
    N_phase = 3,
    population = tab_bins(ages,3)
  )

  # Dose 1
  d1_n <- nimue_target_pop(dose_number = 1, dose_times, prioritisation = rep(1, 3),
                             t = 1, dose_period = 14, d2_prioritise = rep(FALSE, 3))


  d1_s <- target_pop(
    dose = 1,variables = variables,events = events,parameters = parameters,timestep = 1,dt = 1,strategy_matrix_step = rep(1,3)
  )

  d1_s_out <- vapply(d1_s,function(b){b$size()},numeric(1))

  expect_equal(d1_n, d1_s_out)

  # Dose 1 as a function of prioritisation matrix
  d1_pri_n <- nimue_target_pop(dose_number = 1, dose_times, prioritisation = c(0, 1, 0),
                         t = 1, dose_period = 14, d2_prioritise = rep(FALSE, 3))

  d1_pri_s <- target_pop(
    dose = 1,variables = variables,events = events,parameters = parameters,timestep = 1,dt = 1,strategy_matrix_step = c(0,1,0)
  )

  d1_pri_s_out <- vapply(d1_pri_s,function(b){b$size()},numeric(1))

  expect_equal(d1_pri_n, d1_pri_s_out)

  # Dose 2 - none as all d2_prioritise set to FALSE
  d2_n <- nimue_target_pop(dose_number = 2, dose_times, prioritisation = c(1, 1, 1),
                     t = 1, dose_period = 14, d2_prioritise = rep(FALSE, 3))

  d2_s <- safir::target_pop(
    dose = 2,variables = variables,events = events,parameters = parameters,timestep = 1,dt = 1,strategy_matrix_step = rep(1,3),next_dose_priority = rep(0,3)
  )

  d2_s_out <- vapply(d2_s,function(b){b$size()},numeric(1))

  expect_equal(d2_n, d2_s_out)

  # Dose 2 - none as too soon after dose 1
  d2_t_n <- nimue_target_pop(dose_number = 2, dose_times, prioritisation = c(1, 1, 1),
             t = 1, dose_period = 14, d2_prioritise = rep(TRUE, 3))

  d2_t_s <- safir::target_pop(
    dose = 2,variables = variables,events = events,parameters = parameters,timestep = 1,dt = 1,strategy_matrix_step = rep(1,3),next_dose_priority = rep(1,3)
  )

  d2_t_s_out <- vapply(d2_t_s,function(b){b$size()},numeric(1))

  expect_equal(d2_t_n ,d2_t_s_out)

  # Dose 2
  d2_ok_n <- nimue_target_pop(dose_number = 2, dose_times, prioritisation = c(1, 1, 1),
                        t = 15, dose_period = 14, d2_prioritise = rep(TRUE, 3))

  d2_ok_s <- safir::target_pop(
    dose = 2,variables = variables,events = events,parameters = parameters,timestep = 15,dt = 1,strategy_matrix_step = rep(1,3),next_dose_priority = rep(1,3)
  )

  d2_ok_s_out <- vapply(d2_ok_s,function(b){b$size()},numeric(1))

  expect_equal(d2_ok_n, d2_ok_s_out)
})


test_that("target_pop is working in general case", {

  n <- 15
  dose_num <- rep(1,n)

  ages <- rep(1:3,each=5)

  parameters <- list(
    N_age = 3,
    dose_period = c(NaN, 6, 4),
    N_phase = 3,
    population = tab_bins(ages,3),
    correlated = FALSE
  )

  variables <- list()
  variables$discrete_age <- IntegerVariable$new(ages)
  variables <- create_vaccine_variables(variables = variables,parameters = parameters)
  initialize_vaccine_variables(variables = variables,dose_time_init = list(rep(1:3,each=5),rep(-1,n),rep(-1,n)),dose_num_init = dose_num)

  events <- list(scheduled_dose = replicate(n = 3,expr = {TargetedEvent$new(n)}))

  # all phase 1 targets reached
  p1 <- target_pop(
    dose = 1,variables = variables,events = events,parameters = parameters,timestep = 1,dt = 1,strategy_matrix_step = rep(1,3)
  )
  expect_equal(
    vapply(p1,function(b){b$size()},numeric(1)),
    rep(0, 3)
  )

  # at time 8 only groups 1 and 2 are good for phase 2
  p2 <- safir::target_pop(
    dose = 2,variables = variables,events = events,parameters = parameters,timestep = 8,dt = 1,strategy_matrix_step = rep(1,3)
  )
  expect_equal(
    vapply(p2,function(b){b$size()},numeric(1)),
    c(5,5,0)
  )

  # # phase 2 vaccinate group 3 at t = 9, they are prioritized for phase 3
  # variables$dose_time[[2]] <- IntegerVariable$new(c(rep(-1,5), rep(-1,5), rep(9,5)))
  # variables$dose_num <- CategoricalVariable$new(categories = c("0","1","2"),initial_values = c(rep("0",10),rep("1",5)))

  variables$dose_time$queue_update(values = 9, index = 11:15)
  variables$dose_num$queue_update(values = 2, index = 11:15)
  variables$dose_time$.update()
  variables$dose_num$.update()

  # variables$dose_time <- IntegerVariable$new(c(rep(-1,5), rep(-1,5), rep(9,5)))
  # variables$dose_num <- CategoricalVariable$new(categories = c("0","1","2"),initial_values = c(rep("0",10),rep("1",5)))

  # should be up for phase 3 if t =13
  p2_nt <- safir::target_pop(
    dose = 3,variables = variables,events = events,parameters = parameters,timestep = 13,dt = 1,strategy_matrix_step = rep(1,3),next_dose_priority = c(0,0,1)
  )
  expect_equal(
    vapply(p2_nt,function(b){b$size()},numeric(1)),
    c(0,0,5)
  )

  # not if t=12
  p2_t <- safir::target_pop(
    dose = 3,variables = variables,events = events,parameters = parameters,timestep = 12,dt = 1,strategy_matrix_step = rep(1,3),next_dose_priority = c(0,0,1)
  )
  expect_equal(
    vapply(p2_t,function(b){b$size()},numeric(1)),
    c(0,0,0)
  )

  # phase 3, group 3 should all be ready to go regardless of priority if t = 13
  p3_nt <- safir::target_pop(
    dose = 3,variables = variables,events = events,parameters = parameters,timestep = 13,dt = 1,strategy_matrix_step = rep(1,3)
  )
  expect_equal(
    vapply(p3_nt,function(b){b$size()},numeric(1)),
    c(0,0,5)
  )

  #  not if t = 12
  p3_t <- safir::target_pop(
    dose = 3,variables = variables,events = events,parameters = parameters,timestep = 12,dt = 1,strategy_matrix_step = rep(1,3)
  )
  expect_equal(
    vapply(p3_t,function(b){b$size()},numeric(1)),
    c(0,0,0)
  )

})


test_that("assign doses is working for phase 1", {

  parameters <- list(
    population = rep(10,3),
    N_age = 3,
    dose_period = c(NaN, 6, 4),
    N_phase = 3,
    correlated = FALSE
  )

  pop_ages <- rep(1:3,times=parameters$population)

  n <- sum(parameters$population)

  var_local <- list()
  var_local$discrete_age <- IntegerVariable$new(pop_ages)
  var_local <- create_vaccine_variables(variables = var_local,parameters = parameters)
  initialize_vaccine_variables(variables = var_local,dose_time_init = list(rep(-1,n),rep(-1,n),rep(-1,n)),dose_num_init = rep(0,n))

  # won't assign no doses
  events <- list(scheduled_dose = replicate(n = 3,expr = {TargetedEvent$new(n)}))
  targeted <- target_pop(
    dose = 1,variables = var_local,events = events,parameters = parameters,timestep = 14,dt = 1,strategy_matrix_step = rep(1,3)
  )
  assign_doses(
    doses_left = 0,
    events = events,dose = 1,eligible = targeted,parameters = parameters
  )
  expect_equal(
    vapply(X = 1:3,FUN = function(x){events$scheduled_dose[[x]]$get_scheduled()$size()},numeric(1)),
    c(0,0,0)
  )

  # assigns all doses
  events <- list(scheduled_dose = replicate(n = 3,expr = {TargetedEvent$new(n)}))
  targeted <- target_pop(
    dose = 1,variables = var_local,events = events,parameters = parameters,timestep = 14,dt = 1,strategy_matrix_step = rep(1,3)
  )
  assign_doses(
    doses_left = 30,
    events = events,dose = 1,eligible = targeted,parameters = parameters
  )
  sched <- var_local$discrete_age$get_values(events$scheduled_dose[[1]]$get_scheduled())
  expect_equal(
    pop_ages,
    sched
  )

  # assigns partial doses (1)
  events <- list(scheduled_dose = replicate(n = 3,expr = {TargetedEvent$new(n)}))
  targeted <- target_pop(
    dose = 1,variables = var_local,events = events,parameters = parameters,timestep = 14,dt = 1,strategy_matrix_step = rep(1,3)
  )
  assign_doses(
    doses_left = 15,
    events = events,dose = 1,eligible = targeted,parameters = parameters
  )
  sched <- var_local$discrete_age$get_values(events$scheduled_dose[[1]]$get_scheduled())
  expect_equal(
    c(5,5,5),
    as.vector(table(sched))
  )

  # assigns partial doses (2)
  events <- list(scheduled_dose = replicate(n = 3,expr = {TargetedEvent$new(n)}))
  targeted <- target_pop(
    dose = 1,variables = var_local,events = events,parameters = parameters,timestep = 14,dt = 1,strategy_matrix_step = rep(1,3)
  )
  safir::assign_doses(
    doses_left = 16,
    events = events,dose = 1,eligible = targeted,parameters = parameters
  )
  sched <- var_local$discrete_age$get_values(events$scheduled_dose[[1]]$get_scheduled())
  assigned_groups <- as.vector(table(sched))
  expect_equal(
    16, sum(assigned_groups)
  )
  expect_true(
    all(assigned_groups >= 5)
  )
  expect_true(
    all(assigned_groups <= 6)
  )
})


test_that("assign doses is working for phase 2", {

  parameters <- list(
    population = rep(10,3),
    N_age = 3,
    dose_period = c(NaN, 6, 4),
    N_phase = 3,
    correlated = FALSE
  )

  pop_ages <- rep(1:3,times=parameters$population)

  n <- sum(parameters$population)

  var_local <- list()
  var_local$discrete_age <- IntegerVariable$new(pop_ages)
  var_local <- create_vaccine_variables(variables = var_local,parameters = parameters)
  initialize_vaccine_variables(variables = var_local,dose_time_init = list(c(rep(5,n-5),rep(-1,5)),rep(-1,n),rep(-1,n)),dose_num_init = c(rep(1,n-5),rep(0,5)))

  # won't assign doses if not past threshold
  t <- 10
  doses <- 30
  phase <- 2

  events <- list(scheduled_dose = replicate(n = 3,expr = {TargetedEvent$new(n)}))
  targeted <- target_pop(
    dose = phase,variables = var_local,events = events,parameters = parameters,timestep = t,dt = 1,strategy_matrix_step = rep(1,3)
  )
  assign_doses(
    doses_left = doses,
    events = events,dose = phase,eligible = targeted,parameters = parameters
  )

  expect_equal(
    vapply(X = 1:3,FUN = function(x){events$scheduled_dose[[x]]$get_scheduled()$size()},numeric(1)),
    c(0,0,0)
  )

  # will assign doses if past threshold, but only to eligible persons
  t <- 12
  doses <- 30
  phase <- 2

  events <- list(scheduled_dose = replicate(n = 3,expr = {TargetedEvent$new(n)}))
  targeted <- target_pop(
    dose = phase,variables = var_local,events = events,parameters = parameters,timestep = t,dt = 1,strategy_matrix_step = rep(1,3)
  )
  assign_doses(
    doses_left = doses,
    events = events,dose = phase,eligible = targeted,parameters = parameters
  )

  sched_size <- vapply(X = 1:3,FUN = function(x){events$scheduled_dose[[x]]$get_scheduled()$size()},numeric(1))
  sched_age <- var_local$discrete_age$get_values(events$scheduled_dose[[2]]$get_scheduled())
  expect_equal(
    sched_size,
    c(0, 25, 0)
  )
  expect_equal(
    as.vector(table(sched_age)),
    c(10, 10, 5)
  )

  # assign partial priority (2) doses to a certain age group
  t <- 14
  doses <- 8
  phase <- 2

  events <- list(scheduled_dose = replicate(n = 3,expr = {TargetedEvent$new(n)}))
  targeted <- target_pop(
    dose = phase,variables = var_local,events = events,parameters = parameters,timestep = t,dt = 1,strategy_matrix_step = rep(1,3),next_dose_priority = c(1,0,0)
  )
  assign_doses(
    doses_left = doses,
    events = events,dose = phase,eligible = targeted,parameters = parameters
  )

  sched_size <- vapply(X = 1:3,FUN = function(x){events$scheduled_dose[[x]]$get_scheduled()$size()},numeric(1))
  sched_age <- var_local$discrete_age$get_values(events$scheduled_dose[[2]]$get_scheduled())
  expect_equal(
    sched_size,
    c(0, doses, 0)
  )
  expect_equal(
    sched_age,
    rep(1, doses)
  )

})


test_that("allocation does not overallocate", {

  # This test previously would overallocate when using a rmultinom sample as
  # rmultinom is sampling with replacement.

  # You need Multivariate hypergeometric distribution
  # odin:::rmhyper conveniently already written
  parameters <- list(
    population = rep(10,10),
    N_age = 10,
    dose_period = c(NaN, 6, 4),
    N_phase = 3,
    correlated = FALSE
  )

  pop_ages <- rep(1:10,times=parameters$population)

  n <- sum(parameters$population)

  var_local <- list()
  var_local$discrete_age <- IntegerVariable$new(pop_ages)
  var_local <- create_vaccine_variables(variables = var_local,parameters = parameters)
  initialize_vaccine_variables(variables = var_local,dose_time_init = list(rep(-1,n),rep(-1,n),rep(-1,n)),dose_num_init = rep(0,n))


  # assigns partial doses (1)
  events <- list(scheduled_dose = replicate(n = 10,expr = {TargetedEvent$new(n)}))
  targeted <- target_pop(
    dose = 1,variables = var_local,events = events,parameters = parameters,timestep = 14,dt = 1,strategy_matrix_step = rep(1,10)
  )
  # here we end up having to give 9 vaccine doses in the last allocation step
  # which across 10 equal ages rmultinom is likely to allocate more than 1 to a
  # age group
  set.seed(1234L)
  assign_doses(
    doses_left = 99,
    events = events,dose = 1,eligible = targeted,parameters = parameters
  )
  sched <- var_local$discrete_age$get_values(events$scheduled_dose[[1]]$get_scheduled())
  expect_equal(
    c(9,rep(10,9)),
    sort(as.vector(table(sched)))
  )


})
