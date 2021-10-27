test_that("get_time_since_last_dose_or_infection works with both prior infecteds and prior vaccinated", {

  n <- 200

  timestep <- 50
  dt <- 0.5
  dose_num_vec <- sample(x = 0:3,size = n,replace = TRUE)
  inf_num_vec <- sample(x = 0:2,size = n,replace = TRUE)
  who <- vector("list",3)
  when <- vector("list",3)
  N_phase <- 3

  dose_num <- IntegerVariable$new(initial_values = dose_num_vec)
  dose_time <- IntegerVariable$new(initial_values = rep(-1,n))
  for (d in 1:3) {
    who[[d]] <- which(dose_num_vec == d)
    when[[d]] <- sample.int(n = 25,size = length(who[[d]]),replace = TRUE)
    dose_time$queue_update(values = when[[d]],index = who[[d]])
    dose_time$.update()
  }

  inf_num <- IntegerVariable$new(initial_values = inf_num_vec)
  inf_time <- IntegerVariable$new(initial_values = rep(-1, n))
  inf_time$queue_update(values = sample(x = 1:10, size = sum(inf_num_vec == 1), replace = TRUE), index = which(inf_num_vec == 1))
  inf_time$queue_update(values = sample(x = 10:20, size = sum(inf_num_vec == 2), replace = TRUE), index = which(inf_num_vec == 2))
  inf_time$.update()

  vaccinated <- dose_num$get_index_of(set = 0)$not(inplace = TRUE)
  infected <- inf_num$get_index_of(set = 0)$not(inplace = TRUE)

  vaccinated_or_infected <- vaccinated$or(infected)

  # safir
  safir_out <- get_time_since_last_dose_or_infection(timestep = timestep, dt = dt, vaccinated_or_infected = vaccinated_or_infected, dose_time = dose_time, inf_time = inf_time)

  # by hand
  dose_times <- rep(NaN, n)
  for (d in 1:3) {
    dose_times[who[[d]]] <- when[[d]]
  }

  infection_times <- inf_time$get_values()
  infection_times[infection_times < 0] <- NaN

  hand_out <- mapply(FUN = function(dose_tt, inf_tt) {
    if (all(is.nan(dose_tt), is.nan(inf_tt))) {
      return(NULL)
    } else {
      # the latter of the two
      return(pmax(dose_tt, inf_tt, na.rm = TRUE))
    }
  }, dose_tt = dose_times, inf_tt = infection_times, SIMPLIFY = TRUE, USE.NAMES = FALSE)
  hand_out <- unlist(hand_out)
  hand_out <- (timestep - hand_out) * dt

  expect_equal(safir_out, hand_out)

})

test_that("get_time_since_last_dose_or_infection works with both prior vaccinated only", {

  n <- 20

  timestep <- 50
  dt <- 0.5
  dose_num_vec <- sample(x = 0:3,size = n,replace = TRUE)
  inf_num_vec <- rep(0, n)
  who <- vector("list",3)
  when <- vector("list",3)
  N_phase <- 3

  dose_num <- IntegerVariable$new(initial_values = dose_num_vec)
  dose_time <- IntegerVariable$new(initial_values = rep(-1,n))
  for (d in 1:3) {
    who[[d]] <- which(dose_num_vec == d)
    when[[d]] <- sample.int(n = 25,size = length(who[[d]]),replace = TRUE)
    dose_time$queue_update(values = when[[d]],index = who[[d]])
    dose_time$.update()
  }

  inf_num <- IntegerVariable$new(initial_values = inf_num_vec)
  inf_time <- IntegerVariable$new(initial_values = rep(-1, n))

  vaccinated <- dose_num$get_index_of(set = 0)$not(inplace = TRUE)
  infected <- inf_num$get_index_of(set = 0)$not(inplace = TRUE)

  vaccinated_or_infected <- vaccinated$or(infected)

  # safir
  safir_out <- get_time_since_last_dose_or_infection(timestep = timestep, dt = dt, vaccinated_or_infected = vaccinated_or_infected, dose_time = dose_time, inf_time = inf_time)

  # by hand
  hand_out <- rep(NaN,n)
  for (d in 1:3) {
    hand_out[who[[d]]] <- (timestep - when[[d]]) * dt
  }
  hand_out <- hand_out[!is.nan(hand_out)]

  expect_equal(safir_out, hand_out)

})

test_that("get_time_since_last_dose_or_infection works with both prior infected only", {

  n <- 20

  timestep <- 50
  dt <- 0.5
  dose_num_vec <- rep(0, n)
  inf_num_vec <- sample(x = 0:2,size = n,replace = TRUE)
  who <- vector("list",3)
  when <- vector("list",3)
  N_phase <- 3

  dose_num <- IntegerVariable$new(initial_values = dose_num_vec)
  dose_time <- IntegerVariable$new(initial_values = rep(-1,n))

  inf_num <- IntegerVariable$new(initial_values = inf_num_vec)
  inf_time <- IntegerVariable$new(initial_values = rep(-1, n))
  inf_time$queue_update(values = sample(x = 1:10, size = sum(inf_num_vec == 1), replace = TRUE), index = which(inf_num_vec == 1))
  inf_time$queue_update(values = sample(x = 10:20, size = sum(inf_num_vec == 2), replace = TRUE), index = which(inf_num_vec == 2))
  inf_time$.update()

  vaccinated <- dose_num$get_index_of(set = 0)$not(inplace = TRUE)
  infected <- inf_num$get_index_of(set = 0)$not(inplace = TRUE)

  vaccinated_or_infected <- vaccinated$or(infected)

  # safir
  safir_out <- get_time_since_last_dose_or_infection(timestep = timestep, dt = dt, vaccinated_or_infected = vaccinated_or_infected, dose_time = dose_time, inf_time = inf_time)

  # by hand
  infection_times <- inf_time$get_values()
  infection_times <- infection_times[infection_times > 0]
  hand_out <- (timestep - infection_times) * dt

  expect_equal(safir_out, hand_out)

})
