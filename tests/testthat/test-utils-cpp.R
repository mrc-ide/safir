test_that("float equality works", {

  expect_true(compare_floats(1,1.0))
  expect_false(compare_floats(1,0.9999))

})

test_that("cross tab works", {

  a_mar <- 4
  b_mar <- 3

  a <- c(1,2,3,4,1,2,3,4)
  b <- c(1,1,1,2,2,2,3,3)

  expect_true(all(cross_tab_margins(a,b,a_mar,b_mar) == table(a,b)))

})

test_that("cross tab for dose/age works", {

  a <- IntegerVariable$new(0:4)
  b <- IntegerVariable$new(c(1,2,3,1,2))

  comp1 <- cross_tab_doses_age(a$.variable,b$.variable,4,3)
  comp2 <- as.matrix(table(a$get_values(), b$get_values()))

  expect_equal(comp1, matrix(comp2, 5, 3))

  dose_vals <- sample(x = 0:3,size = 100,replace = TRUE)
  age_vals <- sample(x = 1:10,size = 100,replace = TRUE)

  a <- IntegerVariable$new(dose_vals)
  b <- IntegerVariable$new(age_vals)
  comp1 <- cross_tab_doses_age(a$.variable,b$.variable,3,10)
  comp2 <- matrix(table(a$get_values(), b$get_values()), 4, 10)

  expect_equal(comp1, comp2)

  dose_vals <- c(0,0,0,3)
  age_vals <- c(1,2,2,1)
  a <- IntegerVariable$new(dose_vals)
  b <- IntegerVariable$new(age_vals)

  comp1 <- cross_tab_doses_age(a$.variable,b$.variable,3,2)
  comp2 <- matrix(data = c(1, 2, 0, 0, 0, 0, 1, 0),nrow = 4, ncol = 2, byrow = TRUE)
  expect_equal(comp1, comp2)
})

test_that("cross tab for compartments/age works", {

  SIR <- c("S", "I", "R")

  ages <- sample.int(n = 6, size = 100, replace = TRUE)
  compartments <- sample(x = SIR, size = 100, replace = TRUE)
  sir_to_int <- c("S" = 1, "I" = 2, "R" = 3)

  comp1 <- matrix(data = 0, nrow = 6, ncol = 3)
  for (i in 1:100) {
    comp1[ages[i], sir_to_int[[compartments[i]]]] <- comp1[ages[i], sir_to_int[[compartments[i]]]] + 1
  }

  expect_equal(colSums(comp1), as.vector(table(compartments)[SIR]))
  expect_equal(rowSums(comp1), as.vector(table(ages)[as.character(1:6)]))

  comp_variable <- CategoricalVariable$new(categories = SIR, initial_values = compartments)
  age_variable <- IntegerVariable$new(initial_values = ages)

  comp2 <- cross_tab_compartments_age(compartments = comp_variable$.variable, age = age_variable$.variable, num_ages = 6, compartment_names = SIR)

  expect_equal(comp1, comp2)

})


test_that("tab_bins works", {

  nbin <- 10
  a <- sample.int(n = nbin,size = 100,replace = TRUE)

  tabR <- tabulate(bin = a,nbins = nbin)
  tabC <- tab_bins(a = a,nbins = nbin)

  expect_identical(tabR,tabC)

})

test_that("get contact matrix works", {

  iso3c <- "GBR"
  pop <- get_population(iso3c)
  contact_mat <- squire::get_mixing_matrix(iso3c = iso3c)

  parameters <- get_parameters(
    population = pop$n,
    contact_matrix_set = contact_mat,
    iso3c = iso3c,
    time_period = 100,
    dt = 1
  )

  expect_equal(
    get_contact_matrix(parameters = parameters),
    get_contact_matrix_cpp(array = parameters$mix_mat_set,i = 0)
  )

})

test_that("c++ matrix ops work", {

  m <- matrix(rexp(9),3,3)
  mm <- matrix(rexp(9),3,3)
  m1 <- 7:9
  m2 <- 11:13

  expect_equal(
    matrix_vec_mult_cpp(m = m,a = m1),
    as.vector(m %*% m1)
  )

  expect_equal(
    matrix_2vec_mult_cpp(m = m,a = m1,b = m2),
    rowSums(m %*% diag(m1) %*% diag(m2))
  )

  expect_equal(
    mult_2matrix_rowsum(a = m,b = mm),
    rowSums(m * mm)
  )
})

test_that("get proportion vaccinated works", {

  age <- individual::IntegerVariable$new(initial_values = sample.int(10,1e4,T))
  vaxx <- individual::Bitset$new(1e4)
  vaxx$insert(sample.int(n = 1e4,size = 4e3,replace = F))

  pr_r <- vapply(X = 1:10,FUN = function(x){
    get_proportion_vaccinated_nimue(variables = list(discrete_age=age,vaccinated=vaxx),age = x)
  }, numeric(1))

  pr_cpp <- vapply(X = 1:10,FUN = function(x){
    get_proportion_vaccinated_nimue_internal(discrete_age = age$.variable,vaccinated = vaxx$.bitset,age = x)
  }, numeric(1))

  expect_equal(pr_r,pr_cpp)
})
