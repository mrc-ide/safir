test_that("test the correct values are returned from get_asymptomatic", {

  asymp <- get_asymptomatic()

  expect_equal(asymp$dur_IAsymp, 2.1)
  expect_equal(asymp$prob_asymp[4], 0.2)
  expect_equal(length(asymp$IAsymp_0), 17)

})

test_that("test the correct values are returned from get_country", {

  country <- get_country(iso3c = "AFG")

  expect_equal(country, "Afghanistan")

})

test_that("test the correct values are returned from get_population", {

  population <- get_population(iso3c = "AFG")

  expect_equal(population$n[1], 5672509)
  expect_equal(population$n[10], 1304736)
  expect_equal(population$n[17], 105925)

})

test_that("test get_parameters returns the correct values from SQUIRE", {

  R0 <- 2
  timestep <- 100
  time_period <- 1000
  tt_contact_matrix <- 0

  pop <- get_population("AFG")
  pop$n <- as.integer(pop$n / 1000)

  psq <- get_parameters(
    iso3c = "AFG",
    population = pop$n,
    R0 = R0,
    time_period = time_period,
    tt_contact_matrix = tt_contact_matrix,
  )

  expect_equal(psq$dur_E, 4.6)
  expect_equal(psq$N_age, 17)
  expect_equal(length(psq$S_0), 17)
  expect_equal(length(psq$IAsymp_0), 17)
  expect_equal(psq$dur_IAsymp, 2.1)
  expect_equal(psq$prob_asymp[8], 0.2)
  expect_equal(length(psq$IRec1_0), 17)
  expect_equal(psq$time_period, 1000)

})

test_that("test that if max-age is set as NULL get_parameters returns the
          correct values from SQUIRE", {

  R0 <- 2
  timestep <- 100
  time_period <- 1000
  tt_contact_matrix <- 0

  pop <- get_population("AFG")
  pop$n <- as.integer(pop$n / 1000)

  psq <- get_parameters(
    iso3c = "AFG",
    population = pop$n,
    R0 = R0,
    time_period = time_period,
    tt_contact_matrix = tt_contact_matrix,
    max_age = NULL
  )

  expect_equal(psq$dur_E, 4.6)
  expect_equal(psq$N_age, 17)
  expect_equal(length(psq$S_0), 17)
  expect_equal(length(psq$IAsymp_0), 17)
  expect_equal(psq$dur_IAsymp, 2.1)
  expect_equal(psq$prob_asymp[8], 0.2)
  expect_equal(length(psq$IRec1_0), 17)
  expect_equal(psq$time_period, 1000)

})


test_that("test interpolation of pars", {


y <- list(matrix(1:4, 2,2), matrix(5:8,2,2))
x <- c(0,10)
end <- 20

want <- interp_input_par(x, y)


y <- c(1,3)
x <- c(0,10)
end <- 20

want <- interp_input_par(x, y)

})
