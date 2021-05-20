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

# test_that("r_exp works", {
#
#   ret <- r_exp(1, 0.1)
#   expect_true(ret >= 0.0)
#   expect_true(ret <= 1.0)
#
# })

test_that("remove_non_numerics removes characters and characters of arrays", {
  actual <- remove_non_numerics(list(
    a = array(c('1', '2'), dim=c(1, 2)),
    b = c('1', '2'),
    c = c(1, 2)
  ))
  expect_equal(actual, list(c = c(1, 2)))
})

# test_that("bernoulli_multi_p works with empty p", {
#   expect_equal(bernoulli_multi_p(NULL), logical(0))
# })
#
# test_that("bernoulli_multi_p works with two p", {
#   mockery::stub(bernoulli_multi_p, 'runif', mockery::mock(c(.1, .7)))
#   expect_equal(bernoulli_multi_p(c(.2, .6)), c(TRUE, FALSE))
# })
#
# test_that("bernoulli_multi_p works on the boundary", {
#   mockery::stub(bernoulli_multi_p, 'runif', mockery::mock(c(.1, .7)))
#   expect_equal(bernoulli_multi_p(c(.2, .7)), c(TRUE, FALSE))
# })
#
