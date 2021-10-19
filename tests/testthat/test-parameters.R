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
    dt = 1
  )

  expect_equal(attr(psq, "type"), "safir_squire")
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
    max_age = NULL,
    dt = 1
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
expect_equal(dim(want)[3], 11)

y <- c(1,3)
x <- c(0,10)
end <- 20

want <- interp_input_par(x, y)
expect_equal(length(want), 11)

})


test_that("test nimue parameters", {

  iso3c <- "AFG"
  pop <- safir::get_population(iso3c)
  pop$n <- as.integer(pop$n / 1000)
  contact_mat <- squire::get_mixing_matrix(iso3c = iso3c)

  time_period <- 1000
  R0 <- 2
  dt <- 1

  vaccine_coverage_mat <- strategy_matrix(strategy = "Elderly",max_coverage = 0.2)
  tt_vaccine <- c(0, 10:100)
  max_vaccine <- c(0, seq(1e3, 5e4, length.out = length(tt_vaccine)-1))

  parameters <- get_parameters_nimue(
    population = pop$n,
    contact_mat = contact_mat,
    time_period = time_period,
    R0 = R0,
    max_vaccine = max_vaccine,
    tt_vaccine = tt_vaccine,
    vaccine_coverage_mat = vaccine_coverage_mat,
    dt = dt
  )

  expect_equal(attr(parameters, "type"), "safir_nimue")
  expect_equal(parameters$dur_E, 4.6)
  expect_equal(parameters$N_age, 17)
  expect_equal(ncol(parameters$S_0), 6)
  expect_equal(ncol(parameters$S_0), 17)
  expect_equal(length(parameters$IAsymp_0), 17)
  expect_equal(parameters$dur_IAsymp, 2.1)
  expect_equal(parameters$prob_asymp[8], 0.2)
  expect_equal(ncol(parameters$IRec1_0), 6)
  expect_equal(ncol(parameters$IRec1_0), 17)
  expect_equal(parameters$time_period, 1000)
  expect_equal(length(parameters$vaccine_set), time_period+1)

})

test_that("test vaccine titre parameters", {

  ab_pars <- get_vaccine_ab_titre_parameters(vaccine = "Pfizer")
  expect_equal(attr(ab_pars, "type"), "ab_titre")

  ab_pars <- get_vaccine_ab_titre_parameters(vaccine = "Moderna")
  expect_equal(attr(ab_pars, "type"), "ab_titre")

  ab_pars <- get_vaccine_ab_titre_parameters(vaccine = "Sinovac")
  expect_equal(attr(ab_pars, "type"), "ab_titre")

  ab_pars <- get_vaccine_ab_titre_parameters(vaccine = "AstraZeneca")
  expect_equal(attr(ab_pars, "type"), "ab_titre")

  expect_error(get_vaccine_ab_titre_parameters(vaccine = "not_a_vaccine"))

})

test_that("test safir parameters for vaccine (no voc)", {

  iso3c <- "AFG"
  pop <- safir::get_population(iso3c)
  pop$n <- as.integer(pop$n / 1000)
  contact_mat <- squire::get_mixing_matrix(iso3c = iso3c)

  time_period <- 1000
  R0 <- 2
  dt <- 1

  vaccine_doses <- 2
  dose_period <- c(NaN, 28)
  vaccine_set <- vaccine_set <- c(0, seq(from = 1e3, to = 1e4, length.out = time_period-1))
  vaccine_set <- floor(vaccine_set)

  parameters <- safir::get_parameters(
    population = pop$n,
    contact_matrix_set = contact_mat,
    iso3c = iso3c,
    R0 = R0,
    time_period = time_period,
    dt = dt
  )

  ab_parameters <- get_vaccine_ab_titre_parameters(vaccine = "Pfizer", max_dose = vaccine_doses,correlated = TRUE)

  next_dose_priority <- matrix(data = 0, nrow = vaccine_doses - 1,ncol = ncol(vaccine_coverage_mat))
  next_dose_priority[1, 15:17] <- 1 # prioritize 3 oldest age groups for next dose

  # test single strategy matrix
  vaccine_coverage_mat <- strategy_matrix(strategy = "Elderly",max_coverage = 0.5)

  parameters <- make_vaccine_parameters(
    safir_parameters = parameters,
    vaccine_ab_parameters = ab_parameters,
    vaccine_set = vaccine_set,
    dose_period = dose_period,
    strategy_matrix = vaccine_coverage_mat,
    next_dose_priority_matrix = next_dose_priority
  )

  expect_equal(attr(parameters, "type"), "safir_vaccine_notype")
  expect_equal(length(parameters$vaccine_coverage_mat), vaccine_doses)
  expect_true(inherits(parameters$vaccine_coverage_mat, "list"))
  expect_equal(parameters$dur_E, 4.6)
  expect_equal(parameters$N_age, 17)
  expect_equal(length(parameters$S_0), 17)
  expect_equal(length(parameters$IAsymp_0), 17)
  expect_equal(parameters$dur_IAsymp, 2.1)
  expect_equal(parameters$prob_asymp[8], 0.2)
  expect_equal(length(parameters$IRec1_0), 17)
  expect_equal(parameters$time_period, 1000)

  # test list of strategies
  vaccine_coverage_list <- list(strategy_matrix(strategy = "Elderly",max_coverage = 0.5), strategy_matrix(strategy = "Working Elderly Children",max_coverage = 0.8))

  parameters <- make_vaccine_parameters(
    safir_parameters = parameters,
    vaccine_ab_parameters = ab_parameters,
    vaccine_set = vaccine_set,
    dose_period = dose_period,
    strategy_matrix = vaccine_coverage_list,
    next_dose_priority_matrix = next_dose_priority
  )

  expect_equal(attr(parameters, "type"), "safir_vaccine_notype")
  expect_equal(length(parameters$vaccine_coverage_mat), vaccine_doses)
  expect_true(inherits(parameters$vaccine_coverage_mat, "list"))
  expect_equal(parameters$vaccine_coverage_mat, vaccine_coverage_list)
  expect_equal(parameters$dur_E, 4.6)
  expect_equal(parameters$N_age, 17)
  expect_equal(length(parameters$S_0), 17)
  expect_equal(length(parameters$IAsymp_0), 17)
  expect_equal(parameters$dur_IAsymp, 2.1)
  expect_equal(parameters$prob_asymp[8], 0.2)
  expect_equal(length(parameters$IRec1_0), 17)
  expect_equal(parameters$time_period, 1000)

  # failure
  expect_error(
    make_vaccine_parameters(
      safir_parameters = parameters,
      vaccine_ab_parameters = ab_parameters,
      vaccine_set = vaccine_set,
      dose_period = dose_period,
      strategy_matrix = c(vaccine_coverage_list, vaccine_coverage_list),
      next_dose_priority_matrix = next_dose_priority
    )
  )

  expect_error(
    make_vaccine_parameters(
      safir_parameters = parameters,
      vaccine_ab_parameters = ab_parameters,
      vaccine_set = vaccine_set,
      dose_period = dose_period,
      strategy_matrix = matrix(0, nrow = 3, ncol = 3),
      next_dose_priority_matrix = next_dose_priority
    )
  )
})

