test_that("test VFR vector gives reasonable output", {

  # dt = 1

  library(squire)

  iso3c <- "GBR"
  pop <- safir::get_population(iso3c)
  contact_mat <- squire::get_mixing_matrix(iso3c = iso3c)

  # Create our simulation parameters
  time_period <- 200
  dt <- 1
  timesteps <- time_period/dt

  parameters <- get_parameters(
    population = pop$n,
    contact_matrix_set = contact_mat,
    iso3c = iso3c,
    R0 = 2,
    time_period = time_period,
    dt = dt
  )

  vfr_time_1 <- 100
  vfr_time_2 <- 120
  vfr_val <- 4

  vfr <- variant_fold_reduction_vector(parameters = parameters, dt = dt, vfr = vfr_val, vfr_time_1 = vfr_time_1, vfr_time_2 = vfr_time_2)

  expect_equal(length(vfr), timesteps)
  expect_equal(vfr[1:vfr_time_1], rep(1, vfr_time_1))
  expect_true(vfr[vfr_time_1 + 1] < vfr[vfr_time_2])
  expect_equal(vfr[vfr_time_2:timesteps], rep(vfr_val, timesteps - vfr_time_2 + 1))

  # dt = 0.2

  # Create our simulation parameters
  time_period <- 200
  dt <- 0.2
  timesteps <- time_period/dt

  parameters <- get_parameters(
    population = pop$n,
    contact_matrix_set = contact_mat,
    iso3c = iso3c,
    R0 = 2,
    time_period = time_period,
    dt = dt
  )

  vfr_time_1 <- 100
  vfr_time_2 <- 120
  vfr_val <- 4

  vfr <- variant_fold_reduction_vector(parameters = parameters, dt = dt, vfr = vfr_val, vfr_time_1 = vfr_time_1, vfr_time_2 = vfr_time_2)

  expect_equal(length(vfr), timesteps)
  expect_equal(vfr[1:vfr_time_1], rep(1, vfr_time_1))
  expect_true(vfr[(vfr_time_1/dt) + 1] < vfr[vfr_time_2/dt])
  expect_equal(vfr[(vfr_time_2/dt):timesteps], rep(vfr_val, timesteps - (vfr_time_2/dt) + 1))

})


test_that("test VFR decay with NAT from vaccination working with dt = 1", {

  # Create our simulation parameters
  time_period <- 800
  dt <- 1
  n <- 10
  timesteps <- time_period/dt

  vacc_parameters <- get_vaccine_ab_titre_parameters(vaccine = "Pfizer")
  parameters <- make_vaccine_parameters(safir_parameters = list(), vaccine_ab_parameters = vacc_parameters, dose_period = NaN, vaccine_set = rep(0, time_period), strategy_matrix = matrix(0, 17, 17), next_dose_priority_matrix = matrix(0, 0, 17))

  parameters$time_period <- time_period
  parameters$population <- n
  parameters$N_phase <- 1

  vfr_time_1 <- 100
  vfr_time_2 <- 120
  vfr_val <- 4

  # compute with gold standard code
  compare_nt <- draw_nt_vfr(parameters = vacc_parameters, n = n, tmax = time_period, vf = vfr_val, vfr_time_1 = vfr_time_1, vfr_time_2 = vfr_time_2)

  nt <- compare_nt$nt
  z1 <- compare_nt$z1

  # compute with safir
  variables <- create_vaccine_variables(variables = list(), parameters = parameters)
  ef_infection <- matrix(0, nrow = timesteps + 1, ncol = n)
  ef_severe <- matrix(0, nrow = timesteps + 1, ncol = n)
  ef_infection_cpp <- matrix(0, nrow = timesteps + 1, ncol = n)
  ef_severe_cpp <- matrix(0, nrow = timesteps + 1, ncol = n)

  vfr_vector <- variant_fold_reduction_vector(parameters = parameters, dt = dt, vfr = vfr_val, vfr_time_1 = vfr_time_1, vfr_time_2 = vfr_time_2)
  parameters <- make_immune_parameters(parameters = parameters, vfr = vfr_vector)

  # time 0 NAT boost from vaccine
  schedule_dose_vaccine(timestep = 0,variables = variables,target = Bitset$new(n)$insert(1:n),dose = 1,parameters = parameters)
  variables$ab_titre$queue_update(values = log(10^z1)) # make sure using the same set of RVs
  variables$ab_titre$.update()
  variables$dose_num$.update()
  variables$dose_time$.update()

  ab_titre <- vaccine_ab_titre_process(parameters = parameters,variables = variables, dt = dt)

  safir_out <- matrix(data = NaN,nrow = timesteps + 1,ncol = n)
  safir_out[1, ] <- variables$ab_titre$get_values()
  ef_infection[1, ] <- vaccine_efficacy_infection(ab_titre = variables$ab_titre$get_values(), parameters = parameters, timestep = 1)
  ef_severe[1, ] <- vaccine_efficacy_severe(ab_titre = variables$ab_titre$get_values(), ef_infection = ef_infection[1, ], parameters = parameters, timestep = 1)
  ef_infection_cpp[1, ] <- vaccine_efficacy_infection_cpp(ab_titre = variables$ab_titre$get_values(), parameters = parameters, timestep = 0)
  ef_severe_cpp[1, ] <- vaccine_efficacy_severe_cpp(ab_titre = variables$ab_titre$get_values(), ef_infection = ef_infection_cpp[1, ], parameters = parameters, timestep = 0)


  for (t in 1:timesteps) {
    ab_titre(timestep = t)
    variables$ab_titre$.update()
    safir_out[t + 1L, ] <- variables$ab_titre$get_values()
    ef_infection[t + 1L, ] <- vaccine_efficacy_infection(ab_titre = variables$ab_titre$get_values(), parameters = parameters, timestep = t)
    ef_severe[t + 1L, ] <- vaccine_efficacy_severe(ab_titre = variables$ab_titre$get_values(), ef_infection = ef_infection[t + 1L, ], parameters = parameters, timestep = t)
    ef_infection_cpp[t + 1L, ] <- vaccine_efficacy_infection_cpp(ab_titre = variables$ab_titre$get_values(), parameters = parameters, timestep = t - 1)
    ef_severe_cpp[t + 1L, ] <- vaccine_efficacy_severe_cpp(ab_titre = variables$ab_titre$get_values(), ef_infection = ef_infection_cpp[t + 1L, ], parameters = parameters, timestep = t - 1)

  }

  # test ln scale NAT same
  expect_equal(safir_out, nt)
  expect_equal(1 - ef_infection, compare_nt$ef_infection)
  expect_equal(1 - ef_severe, compare_nt$ef_severe)
  expect_equal(1 - ef_infection_cpp, compare_nt$ef_infection)
  expect_equal(1 - ef_severe_cpp, compare_nt$ef_severe)

})


# test_that("test VFR decay working with dt = 0.5", {
#
#
#
# })
