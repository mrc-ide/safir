test_that("get contact matrix works", {

  n <- 1e3
  N_age <- 5
  doses <- 3
  num_dose <- rbinom(n = 3,size = n,prob = 0.5)

  variables <- list(
    discrete_age = IntegerVariable$new(initial_values = sample.int(n = N_age,size = n,replace = TRUE)),
    dose_time = replicate(n = doses,expr = {IntegerVariable$new(initial_values = rep(-1,n))})
  )

  for (i in 1:doses) {
    variables$dose_time[[i]]$queue_update(values = pmax(1,rpois(n = num_dose[i],lambda = 5)),index = sample.int(n = n,size = num_dose[i],replace = FALSE))
    variables$dose_time[[i]]$.update()
  }

  for (i in 1:doses) {

    R_out <- sapply(X = 1:N_age,FUN = function(a){
      get_proportion_vaccinated(variables = variables,age = a,dose = i)
    })

    cpp_out <- get_proportion_vaccinated_all_ages_cpp(variables = variables,N_age = N_age,dose = i)

    expect_equal(R_out, cpp_out)

  }

})
