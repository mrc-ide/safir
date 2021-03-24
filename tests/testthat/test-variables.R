test_that("create_variables returns the correct output", {

  pop <- get_population("AFG")

  pop$n <- as.integer(pop$n/10000)
  theages <- create_variables(pop, get_parameters("AFG"))
  expect_length(length(theages$age), 1)

  expect_length(theages$age$initial_values, sum(pop$n))

})


test_that("create_continuous_age_variable creates the right number of ages", {

  pop <- get_population("AFG")
  pop$n <- as.integer(pop$n/10000)
  ages <- create_continuous_age_variable(pop)
  expect_length(ages, sum(pop$n))

})

test_that("test create_continuous_age_variable", {

  pop <- squire::get_population(iso3c = "ATG", simple_SEIR = FALSE)
  pop$n <- as.integer(pop$n/100)
  age_cont <- create_continuous_age_variable(pop)

  expect_equal(length(age_cont), sum(pop$n))

})

test_that("test create_discrete_age_variable", {

  pop <- squire::get_population(iso3c = "ATG", simple_SEIR = FALSE)
  pop$n <- as.integer(pop$n/100)
  ages <- create_continuous_age_variable(pop = pop)
  disc_ages <- create_discrete_age_variable(ages, pop)

  expect_equal(as.numeric(table(disc_ages)), pop$n)
})

test_that("test identify_ages_to_adjust and swap_ages work", {

  # Create our parameters
  pop <- squire::get_population(iso3c = "ATG")

  pop$n <- as.integer(pop$n) / 100

  parameters <- get_parameters(
    iso3c = "ATG", population = pop$n
  )

  cont_age <- create_continuous_age_variable(pop, parameters$max_age)

  discrete_age <- create_discrete_age_variable(cont_age, pop)

  swaps <- identify_ages_to_adjust(discrete_age, parameters)

  actual <- swap_ages(swaps, discrete_age)

  # Check that values agree
  e1 <- parameters$E1_0

  expect_equal(
    tail(actual, 20),
    rep(which(e1 > 0), e1[e1 > 0])
  )

})
