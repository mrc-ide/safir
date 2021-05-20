test_that('get_contact_matrices returns the first matrix', {
  data <- c(rep(0, 9), rep(1, 9)) # Matrix 1 and Matrix 2
  m <- aperm(array(data, dim=c(3, 3, 2)), c(3, 1, 2)) # Organise it how R does
  pars <- list(mix_mat_set=m)
  expect_equal(
    get_contact_matrix(pars),
    array(rep(0, 9), dim=c(3, 3))
  )
})
