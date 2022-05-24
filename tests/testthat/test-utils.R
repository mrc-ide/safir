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

test_that("interpolating Rt works", {
  expect_error(interpolate_rt(dates = c(as.Date(x = "2/1/2020", format = "%m/%d/%Y"), as.Date(x = "1/1/2020", format = "%m/%d/%Y")),rt = c(1, 2)))
  expect_error(interpolate_rt(dates = c(as.Date(x = "2/1/2020", format = "%m/%d/%Y"), as.Date(x = "1/1/2020", format = "%m/%d/%Y")),rt = c(1, 2), max_date = as.Date(x = "6/1/2020", format = "%m/%d/%Y")))
  expect_error(interpolate_rt(dates = c(as.Date(x = "2/1/2020", format = "%m/%d/%Y"), as.Date(x = "4/1/2020", format = "%m/%d/%Y")),rt = c(1)))
  expect_error(interpolate_rt(dates = c(as.Date(x = "2/1/2020", format = "%m/%d/%Y"), as.Date(x = "4/1/2020", format = "%m/%d/%Y")),rt = c(1), max_date = as.Date(x = "6/1/2020", format = "%m/%d/%Y")))
  expect_error(interpolate_rt(dates = as.Date(x = "2/1/2020", format = "%m/%d/%Y"),rt = c(1, 2)))
  expect_error(interpolate_rt(dates = as.Date(x = "2/1/2020", format = "%m/%d/%Y"),rt = c(1, 2), max_date = as.Date(x = "6/1/2020", format = "%m/%d/%Y")))
  expect_error(interpolate_rt(dates = "2/1/2020",rt = c(1, 2)))
  expect_error(interpolate_rt(dates = c("2/1/2020", "4/1/2020"),rt = c(1, 2)))
  expect_error(interpolate_rt(dates = as.Date(c("2/1/2020", "4/1/2020"), format = "%m/%d/%Y") ,rt = c(NaN, 2)))
  expect_error(interpolate_rt(dates = as.Date(c("2/1/2020", "4/1/2020"), format = "%m/%d/%Y") ,rt = c(1, -5)))
  expect_error(interpolate_rt(dates = as.Date(c("2/1/2020", "4/1/2020"), format = "%m/%d/%Y") ,rt = c(1, NA)))
  expect_error(interpolate_rt(dates = as.Date(c("2/1/2020", "4/1/2020"), format = "%m/%d/%Y") ,rt = c("5", 2)))

  interp <- interpolate_rt(dates = as.Date(c("2/1/2020", "2/15/2020"), format = "%m/%d/%Y") ,rt = c(1, 2))
  expect_true(all(vapply(interp, length, integer(1)) == 15))
  expect_true(interp$Rt[1] == 1)
  expect_true(interp$Rt[15] == 2)
  expect_equal(interp$Rt_tt, 1:15)
})

test_that("get_size_bset works", {

  max <- rpois(n = 10, lambda = 100)
  fill <- rbinom(n = length(max), size = max, prob = 0.5)
  bsets <- lapply(X = 1:length(max), FUN = function(i) {
    Bitset$new(size = max[i])$insert(1:fill[i])
  })
  expect_equal(get_size_bset(bsets = bsets), fill)
})
