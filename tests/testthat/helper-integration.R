mock_category <- function(...) {
  v <- individual::CategoricalVariable$new(...)
  list(
    get_index_of = v$get_index_of,
    get_size_of = v$get_size_of,
    queue_update = mockery::mock()
  )
}

mock_double <- function(...) {
  v <- individual::DoubleVariable$new(...)
  list(
    get_values = v$get_values,
    queue_update = mockery::mock()
  )
}

mock_integer <- function(...) {
  v <- individual::IntegerVariable$new(...)
  list(
    get_values = v$get_values,
    queue_update = mockery::mock()
  )
}

mock_render <- function(...) {
  v <- individual::Render$new(...)
  list(
    render = mockery::mock()
  )
}

mock_event <- function(event) {
  list(
    get_scheduled = function(...) event$get_scheduled(...),
    schedule = mockery::mock(),
    clear_schedule = mockery::mock()
  )
}

expect_bitset_update <- function(mock, value, index, call = 1) {
  expect_equal(mockery::mock_args(mock)[[call]][[1]], value)
  if("Bitset" %in% is(mockery::mock_args(mock)[[call]][[2]])){
    expect_equal(mockery::mock_args(mock)[[call]][[2]]$to_vector(), index)
  } else {
    expect_equal(mockery::mock_args(mock)[[call]][[2]], index)
  }
}
