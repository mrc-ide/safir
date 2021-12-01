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



test_that("test VFR decay working with dt = 1", {



})


test_that("test VFR decay working with dt = 0.5", {



})
