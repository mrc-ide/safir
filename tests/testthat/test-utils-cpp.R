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

test_that("get proportion vaccinated works", {
  age <- individual::IntegerVariable$new(initial_values = sample.int(10,1e4,T))
  vaxx <- individual::Bitset$new(1e4)
  vaxx$insert(sample.int(n = 1e4,size = 4e3,replace = F))

  pr_r <- sapply(X = 1:10,FUN = function(x){
    get_proportion_vaccinated_nimue(variables = list(discrete_age=age,vaccinated=vaxx),age = x)
  })

  pr_cpp <- sapply(X = 1:10,FUN = function(x){
    get_proportion_vaccinated_nimue_internal(discrete_age = age$.variable,vaccinated = vaxx$.bitset,age = x)
  })

  expect_equal(pr_r,pr_cpp)
})
