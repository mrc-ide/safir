test_that("get_proportion_vaccinated_all_ages_cpp works", {

  library(individual)

  n <- 1e3
  N_age <- 5
  doses <- 3
  num_dose <- rbinom(n = 3,size = n,prob = 0.5)
  ages <- sample.int(n = N_age,size = n,replace = TRUE)
  ages_size <- tab_bins(a = ages,nbins = N_age)

  variables <- list(
    discrete_age = IntegerVariable$new(initial_values = ages),
    dose_time = replicate(n = doses,expr = {IntegerVariable$new(initial_values = rep(-1,n))})
  )

  events <- list(
    scheduled_dose = replicate(n = doses,expr = {TargetedEvent$new(population_size = n)})
  )

  parameters <- list(N_age = N_age)

  for (i in 1:doses) {
    variables$dose_time[[i]]$queue_update(values = pmax(1,rpois(n = num_dose[i],lambda = 5)),index = sample.int(n = n,size = num_dose[i],replace = FALSE))
    variables$dose_time[[i]]$.update()
  }

  for (i in 1:doses) {

    R_out <- get_current_coverage(variables = variables,events = events,dose = i,parameters = parameters)
    R_out <- vapply(R_out, function(a){a$size()}, numeric(1))
    R_out <- R_out / ages_size

    cpp_out <- get_proportion_vaccinated_all_ages_cpp(variables = variables,N_age = N_age,dose = i)

    expect_equal(R_out, cpp_out)

  }

})
