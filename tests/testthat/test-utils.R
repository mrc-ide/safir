context("utils")

test_that("null-or-value works", {
  expect_equal(1 %||% NULL, 1)
  expect_equal(1 %||% 2, 1)
  expect_equal(NULL %||% NULL, NULL)
  expect_equal(NULL %||% 2, 2)
})

test_that("vcapply works", {

  expect_length(vcapply(1:5, function(x) {letters[x]}), 5)
  expect_equal(vcapply(1:5, function(x) {letters[x]}), letters[1:5])

})

test_that("remove_non_numerics removes characters and characters of arrays", {
  actual <- remove_non_numerics(list(
    a = array(c('1', '2'), dim=c(1, 2)),
    b = c('1', '2'),
    c = c(1, 2)
  ))
  expect_equal(actual, list(c = c(1, 2)))
})

test_that("cross tab works", {

  a_mar <- 4
  b_mar <- 3

  a <- c(1,2,3,4,1,2,3,4)
  b <- c(1,1,1,2,2,2,3,3)

  expect_true(all(cross_tab_margins(a,b,a_mar,b_mar) == table(a,b)))

})


