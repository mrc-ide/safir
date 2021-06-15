test_that("coverage and get_proportion_vaccinated are giving the same results", {
  dose_times <- list(matrix(c(1, 2, NA, 2, 3, NA), nrow = 3),
                     matrix(c(NA, NA, NA, 2, 3, 4), nrow = 3),
                     matrix(c(1, 2, 2, NA, NA, NA), nrow = 3))

  dose_1 <- unlist(lapply(dose_times,function(x){x[,1]}))
  dose_1[which(is.na(dose_1))] <- -1

  dose_2 <- unlist(lapply(dose_times,function(x){x[,2]}))
  dose_2[which(is.na(dose_2))] <- -1

  variables <- list()
  variables$dose_time <- list()
  variables$dose_type <- list()
  variables$dose_time[[1]] <- IntegerVariable$new(dose_1)
  variables$dose_time[[2]] <- IntegerVariable$new(dose_2)
  variables$dose_type[[1]] <- CategoricalVariable$new(categories = c("a"),initial_values = rep("a",length(dose_1)))
  variables$dose_type[[2]] <- CategoricalVariable$new(categories = c("a"),initial_values = rep("a",length(dose_2)))
  variables$discrete_age <- IntegerVariable$new(rep(1:length(dose_times),times=sapply(dose_times,nrow)))

  cov_safir <- sapply(X = 1:length(dose_times),FUN = function(a){
    get_proportion_vaccinated(variables = variables,age = a,dose = 1)
  })

  cov_safir_2 <- sapply(X = 1:length(dose_times),FUN = function(a){
    get_proportion_vaccinated(variables = variables,age = a,dose = 2)
  })

  cov_safir_type <- sapply(X = 1:length(dose_times),FUN = function(a){
    get_proportion_vaccinated_type(variables = variables,age = a,dose = 1,type = "a")
  })

  cov_safir_type_2 <- sapply(X = 1:length(dose_times),FUN = function(a){
    get_proportion_vaccinated_type(variables = variables,age = a,dose = 2,type = "a")
  })

  expect_equal(
    cov_safir,
    nimue:::coverage(dose_times, 1)
  )
  expect_equal(
    cov_safir_2,
    nimue:::coverage(dose_times, 2)
  )
  expect_equal(
    cov_safir_type,
    nimue:::coverage(dose_times, 1)
  )
  expect_equal(
    cov_safir_type_2,
    nimue:::coverage(dose_times, 2)
  )
})


test_that('eligable_for_second and eligible_for_dose_vaccine give equivalent results with same input', {
  dose_times <- list(matrix(c(1, 2, NA, NA, 3, NA), nrow = 3),
                     matrix(c(NA, NA, NA, NA, 3, 4), nrow = 3),
                     matrix(c(1, 2, 2, NA, NA, NA), nrow = 3))

  dose_1 <- unlist(lapply(dose_times,function(x){x[,1]}))
  dose_1[which(is.na(dose_1))] <- -1

  dose_2 <- unlist(lapply(dose_times,function(x){x[,2]}))
  dose_2[which(is.na(dose_2))] <- -1

  variables <- list()
  variables$dose_time <- list()
  variables$dose_time[[1]] <- IntegerVariable$new(dose_1)
  variables$dose_time[[2]] <- IntegerVariable$new(dose_2)

  t <- 1
  dt <- 1

  parameters <- list(
    dose_period = c(NaN, 14),
    N_age = 3
  )

  eligible <- eligible_for_dose_vaccine(dose = 2,parameters = parameters,variables = variables,t = t,dt = dt)

  expect_equal(
    which(unlist(nimue:::eligable_for_second(dose_times, t, parameters$dose_period[2]))),
    eligible$to_vector()
  )

  t <- 1
  parameters$dose_period[2] <- 0
  expect_equal(
    eligible_for_dose_vaccine(dose = 2,parameters = parameters,variables = variables,t = t,dt = dt)$to_vector(),
    which(unlist(nimue:::eligable_for_second(dose_times, t, 0)))
  )

  t <- 200
  parameters$dose_period[2] <- 14
  expect_equal(
    eligible_for_dose_vaccine(dose = 2,parameters = parameters,variables = variables,t = t,dt = dt)$to_vector(),
    which(unlist(nimue:::eligable_for_second(dose_times, t, 14)))
  )
})

test_that('eligable_for_second and age_group_eligible_for_dose_vaccine give equivalent results with same input', {
  dose_times <- list(matrix(c(1, 2, NA, NA, 3, NA), nrow = 3),
                     matrix(c(NA, NA, NA, NA, 3, 4), nrow = 3),
                     matrix(c(1, 2, 2, NA, NA, NA), nrow = 3))

  dose_1 <- unlist(lapply(dose_times,function(x){x[,1]}))
  dose_1[which(is.na(dose_1))] <- -1

  dose_2 <- unlist(lapply(dose_times,function(x){x[,2]}))
  dose_2[which(is.na(dose_2))] <- -1

  variables <- list()
  variables$discrete_age <- IntegerVariable$new(rep(1:3,each=3))
  variables$dose_time <- list()
  variables$dose_time[[1]] <- IntegerVariable$new(dose_1)
  variables$dose_time[[2]] <- IntegerVariable$new(dose_2)

  t <- 200
  dt <- 1

  parameters <- list(
    dose_period = c(NaN, 14),
    N_age = 3
  )

  eligible <- age_group_eligible_for_dose_vaccine(dose = 2,parameters = parameters,variables = variables,t = t,dt = dt)

  expect_equal(
    sapply(nimue:::eligable_for_second(dose_times, t, parameters$dose_period[2]),sum),
    eligible
  )
})


test_that('prioritization steps are working', {

  n <- 1e4
  ages_size <- distribute(n = n,p = 17)
  ages <- rep(x = 1:17, times = ages_size)

  strat <- strategy_matrix(strategy = "Elderly")

  parameters <- list(
    N_age = 17,
    N_prioritisation_steps = nrow(strat),
    vaccine_coverage_mat = strat
  )

  for (i in 1:nrow(strat)) {

    variables <- list(discrete_age = IntegerVariable$new(initial_values = ages), dose_time = NULL)
    variables$dose_time[[1]] <- IntegerVariable$new(initial_values = rep(-1,n))

    ages_to_vaxx <- which(strat[i, ] > 0)
    perc_to_vaxx <- strat[i, ages_to_vaxx]
    for (a in seq_along(ages_to_vaxx)) {
      bset_a <- variables$discrete_age$get_index_of(set = ages_to_vaxx[a])
      bset_a <- filter_bitset(bitset = bset_a,other = 1:floor(bset_a$size() * (perc_to_vaxx[a]+0.05)))
      variables$dose_time[[1]]$queue_update(values = 1,index = bset_a)
      variables$dose_time[[1]]$.update()
    }

    step <- get_current_prioritization_step(variables = variables,parameters = parameters,dose = 1)
    if (i < nrow(strat)) {
      expect_equal(
        step, i + 1
      )
    } else {
      expect_equal(
        step, i
      )
    }
  }

})


test_that("target_pop are giving the same results between safir and nimue", {
  dose_times <- list(matrix(c(1, 2, NA, NA, 3, NA), nrow = 3),
                     matrix(c(NA, NA, NA, NA, 3, 4), nrow = 3),
                     matrix(c(1, 2, 2, NA, NA, NA), nrow = 3))

  dose_1 <- unlist(lapply(dose_times,function(x){x[,1]}))
  dose_1[which(is.na(dose_1))] <- -1

  dose_2 <- unlist(lapply(dose_times,function(x){x[,2]}))
  dose_2[which(is.na(dose_2))] <- -1

  variables <- list()
  variables$dose_time <- list()
  variables$dose_type <- list()
  variables$dose_time[[1]] <- IntegerVariable$new(dose_1)
  variables$dose_time[[2]] <- IntegerVariable$new(dose_2)
  variables$discrete_age <- IntegerVariable$new(rep(1:length(dose_times),times=sapply(dose_times,nrow)))

  parameters <- list(
    N_age = 3,
    dose_period = c(NaN, 14, NaN),
    N_phase = 3
  )

  # Dose 1
  d1_n <- nimue:::target_pop(dose_number = 1, dose_times, prioritisation = rep(1, 3),
                     t = 1, dose_period = 14, d2_prioritise = rep(FALSE, 3))

  variables$phase <- 1
  d1_s <- safir::target_pop(
    phase = 1,variables = variables,parameters = parameters,t = 1,dt = 1,prioritisation = rep(1,3)
  )

  expect_equal(d1_n,d1_s)

  # Dose 1 as a function of prioritisation matrix
  d1_pri_n <- nimue:::target_pop(dose_number = 1, dose_times, prioritisation = c(0, 1, 0),
                         t = 1, dose_period = 14, d2_prioritise = rep(FALSE, 3))

  d1_pri_s <- safir::target_pop(
    phase = 1,variables = variables,parameters = parameters,t = 1,dt = 1,prioritisation = c(0,1,0)
  )

  expect_equal(d1_pri_n,d1_pri_s)

  # Dose 2 - none as all d2_prioritise set to FALSE
  d2_n <- nimue:::target_pop(dose_number = 2, dose_times, prioritisation = c(1, 1, 1),
                     t = 1, dose_period = 14, d2_prioritise = rep(FALSE, 3))

  variables$phase <- 1
  d2_s <- safir::target_pop(
    phase = 2,variables = variables,parameters = parameters,t = 1,dt = 1,prioritisation = rep(1,3),vaxx_priority = rep(0, 3)
  )

  expect_equal(d2_n,d2_s)

  # Dose 2 - none as too soon after dose 1
  d2_t_n <- nimue:::target_pop(dose_number = 2, dose_times, prioritisation = c(1, 1, 1),
             t = 1, dose_period = 14, d2_prioritise = rep(TRUE, 3))

  variables$phase <- 1
  d2_t_s <- safir::target_pop(
    phase = 2,variables = variables,parameters = parameters,t = 1,dt = 1,prioritisation = rep(1,3),vaxx_priority = rep(1, 3)
  )

  expect_equal(d2_t_n,d2_t_s)

  # Dose 2
  d2_ok_n <- nimue:::target_pop(dose_number = 2, dose_times, prioritisation = c(1, 1, 1),
                        t = 15, dose_period = 14, d2_prioritise = rep(TRUE, 3))

  variables$phase <- 1
  d2_ok_s <- safir::target_pop(
    phase = 2,variables = variables,parameters = parameters,t = 15,dt = 1,prioritisation = rep(1,3),vaxx_priority = rep(1, 3)
  )

  expect_equal(d2_ok_n,d2_ok_s)
})


test_that("target_pop working in general case", {

  n <- 15
  variables <- list(
    discrete_age = IntegerVariable$new(rep(1:3,each=5)),
    dose_time = NULL
  )
  variables$dose_time[[1]] <- IntegerVariable$new(rep(1:3,each=5))
  variables$dose_time[[2]] <- IntegerVariable$new(rep(-1,n))
  variables$dose_time[[3]] <- IntegerVariable$new(rep(-1,n))

  parameters <- list(
    N_age = 3,
    dose_period = c(NaN, 6, 4),
    N_phase = 3
  )

  # all phase 1 targets reached
  variables$phase <- 1
  expect_equal(
    safir::target_pop(
      phase = 1,variables = variables,parameters = parameters,t = 1,dt = 1,prioritisation = rep(1,3)
    ),
    rep(0, 3)
  )

  # at time 8 only groups 1 and 2 are good for phase 2
  variables$phase <- 2
  expect_equal(
    safir::target_pop(
      phase = 2,variables = variables,parameters = parameters,t = 8,dt = 1,prioritisation = rep(1,3)
    ),
    c(5,5,0)
  )

  # phase 2 vaccinate group 3 at t = 9, they are prioritized for phase 3
  variables$dose_time[[2]] <- IntegerVariable$new(c(rep(-1,5), rep(-1,5), rep(9,5)))
  variables$phase <- 2
  # should be up for phase 3 if t =13
  expect_equal(
    safir::target_pop(
      phase = 3,variables = variables,parameters = parameters,t = 13,dt = 1,prioritisation = rep(1,3),vaxx_priority = c(0,0,1)
    ),
    c(0,0,5)
  )
  # not if t=12
  expect_equal(
    safir::target_pop(
      phase = 3,variables = variables,parameters = parameters,t = 12,dt = 1,prioritisation = rep(1,3),vaxx_priority = c(0,0,1)
    ),
    c(0,0,0)
  )

  # phase 3, group 3 should all be ready to go regardless of priority if t = 13
  variables$phase <- 3
  expect_equal(
    safir::target_pop(
      phase = 3,variables = variables,parameters = parameters,t = 13,dt = 1,prioritisation = rep(1,3)
    ),
    c(0,0,5)
  )
  #  not if t = 12
  expect_equal(
    safir::target_pop(
      phase = 3,variables = variables,parameters = parameters,t = 12,dt = 1,prioritisation = rep(1,3)
    ),
    c(0,0,0)
  )

})
