test_that("cross tab works", {

  a_mar <- 4
  b_mar <- 3

  a <- c(1,2,3,4,1,2,3,4)
  b <- c(1,1,1,2,2,2,3,3)

  expect_true(all(cross_tab_margins(a,b,a_mar,b_mar) == table(a,b)))

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
