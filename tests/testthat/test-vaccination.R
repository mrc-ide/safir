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
  variables$dose_time[[1]] <- IntegerVariable$new(dose_1)
  variables$dose_time[[2]] <- IntegerVariable$new(dose_2)
  variables$discrete_age <- IntegerVariable$new(rep(1:length(dose_times),times=sapply(dose_times,nrow)))

  cov_safir <- sapply(X = 1:length(dose_times),FUN = function(a){
    get_proportion_vaccinated(variables = variables,age = a,dose = 1)
  })

  cov_safir_2 <- sapply(X = 1:length(dose_times),FUN = function(a){
    get_proportion_vaccinated(variables = variables,age = a,dose = 2)
  })

  expect_equal(
    cov_safir,
    nimue:::coverage(dose_times, 1)
  )
  expect_equal(
    cov_safir_2,
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

  t <- 200
  dose_period <- 14
  dt <- 1

  eligible <- eligible_for_dose_vaccine(dose_number = 2,dose_period = dose_period,variables = variables,t = t,dt = dt)

  expect_equal(
    which(unlist(nimue:::eligable_for_second(dose_times, t, dose_period))),
    eligible$to_vector()
  )
})
