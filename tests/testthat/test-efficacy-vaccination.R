test_that("get_time_since_last_dose works", {

  n <- 100

  timestep <- 50
  dt <- 0.5
  dose_num_vec <- sample(x = 0:3,size = n,replace = TRUE)
  who <- vector("list",3)
  when <- vector("list",3)

  N_phase <- 3

  dose_num <- IntegerVariable$new(initial_values = dose_num_vec)
  dose_time <- replicate(n = 3,expr = {IntegerVariable$new(initial_values = rep(-1,n))})
  for (d in 1:3) {
    who[[d]] <- which(dose_num_vec == d)
    when[[d]] <- sample.int(n = 25,size = length(who[[d]]),replace = TRUE)
    dose_time[[d]]$queue_update(values = when[[d]],index = who[[d]])
    dose_time[[d]]$.update()
  }

  vaccinated <- dose_num$get_index_of(set = 0)
  vaccinated <- vaccinated$not()

  # safir
  safir_out <- get_time_since_last_dose(timestep = timestep,dt = dt,vaccinated = vaccinated,dose_num = dose_num,dose_time = dose_time,N_phase = N_phase)

  # by hand
  hand_out <- rep(NaN,n)
  for (d in 1:3) {
    hand_out[who[[d]]] <- (timestep - when[[d]]) * dt
  }
  hand_out <- hand_out[!is.nan(hand_out)]

  expect_equal(safir_out, hand_out)

})
