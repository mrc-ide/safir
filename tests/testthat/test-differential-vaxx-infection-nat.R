test_that("differential decay for infection and vaccine derived NATs works, dt = 1", {

  # setup pop
  dt <- 1
  pop <- 1e2
  pop_assign <- sample(x = c("v", "i" , "vi", "na"), size = pop, replace = TRUE)

  who_any_v <- which(pop_assign %in% c("v", "vi"))
  who_any_i <- which(pop_assign %in% c("i", "vi"))

  nat_v <- rep(-Inf, pop)
  nat_v[who_any_v] <- rexp(length(who_any_v))

  nat_i <- rep(-Inf, pop)
  nat_i[who_any_i] <- rexp(length(who_any_i))

  dose_num0 <- rep(0, pop)
  dose_num0[who_any_v] <- sample(x = 1:2, size = length(who_any_v), replace = TRUE)

  dose_time0 <- rep(-1L, pop)
  dose_time0[who_any_v] <- ceiling(runif(n = length(who_any_v), min = 0, max = 19)/dt)

  inf_num0 <- rep(0, pop)
  inf_num0[who_any_i] <- sample(x = 1:5, size = length(who_any_i), replace = TRUE)

  inf_time0 <- rep(-1L, pop)
  inf_time0[who_any_i] <- ceiling(runif(n = length(who_any_i), min = 0, max = 19)/dt)

  dr_vec_doses <- cbind(runif(n = 40, min = -0.25, max = -0.1), runif(n = 40, min = -0.25, max = -0.1))
  dr_vec_inf <- runif(n = 40, min = -0.25, max = -0.1)

  day <- 20
  timestep <- day/dt

  # manual update: vaccine NATs
  nat_v_update <- nat_v

  who_any_v_time_since_dose <- (timestep - dose_time0[who_any_v]) * dt
  who_any_v_time_since_dose <- ceiling(who_any_v_time_since_dose)

  sub_mat <- cbind(who_any_v_time_since_dose, dose_num0[who_any_v])

  nat_v_update[who_any_v] <- nat_v[who_any_v] + (dr_vec_doses[sub_mat] * dt)

  # manual update: infection NATs
  nat_i_update <- nat_i

  who_any_i_time_since_dose <- (timestep - inf_time0[who_any_i]) * dt
  who_any_i_time_since_dose <- ceiling(who_any_i_time_since_dose)

  nat_i_update[who_any_i] <- nat_i[who_any_i] + (dr_vec_inf[who_any_i_time_since_dose] * dt)

  # safir update
  variables <- list(
    ab_titre = DoubleVariable$new(initial_values = nat_v),
    ab_titre_inf = DoubleVariable$new(initial_values = nat_i),
    dose_num = IntegerVariable$new(initial_values = dose_num0),
    dose_time = IntegerVariable$new(initial_values = dose_time0),
    inf_num = IntegerVariable$new(initial_values = inf_num0),
    inf_time = IntegerVariable$new(initial_values = inf_time0)
  )

  parameters <- list(
    dr_vec_doses = dr_vec_doses,
    dr_vec_inf = dr_vec_inf
  )

  proc <- independent_ab_titre_process(parameters = parameters, variables = variables, dt = dt)
  proc(timestep = timestep)
  for (v in variables) {
    v$.update()
  }

  nat_v_update_safir <- variables$ab_titre$get_values()
  nat_i_update_safir <- variables$ab_titre_inf$get_values()

  expect_equal(nat_v_update, nat_v_update_safir)
  expect_equal(nat_i_update, nat_i_update_safir)

})


test_that("differential decay for infection and vaccine derived NATs works, no vi, dt = 0.5", {

  # setup pop
  dt <- 0.5
  pop <- 1e2
  pop_assign <- sample(x = c("v", "i" , "vi", "na"), size = pop, replace = TRUE)

  who_any_v <- which(pop_assign %in% c("v", "vi"))
  who_any_i <- which(pop_assign %in% c("i", "vi"))

  nat_v <- rep(-Inf, pop)
  nat_v[who_any_v] <- rexp(length(who_any_v))

  nat_i <- rep(-Inf, pop)
  nat_i[who_any_i] <- rexp(length(who_any_i))

  dose_num0 <- rep(0, pop)
  dose_num0[who_any_v] <- sample(x = 1:2, size = length(who_any_v), replace = TRUE)

  dose_time0 <- rep(-1L, pop)
  dose_time0[who_any_v] <- ceiling(runif(n = length(who_any_v), min = 0, max = 19)/dt)

  inf_num0 <- rep(0, pop)
  inf_num0[who_any_i] <- sample(x = 1:5, size = length(who_any_i), replace = TRUE)

  inf_time0 <- rep(-1L, pop)
  inf_time0[who_any_i] <- ceiling(runif(n = length(who_any_i), min = 0, max = 19)/dt)

  dr_vec_doses <- cbind(runif(n = 40, min = -0.25, max = -0.1), runif(n = 40, min = -0.25, max = -0.1))
  dr_vec_inf <- runif(n = 40, min = -0.25, max = -0.1)

  day <- 20
  timestep <- day/dt

  # manual update: vaccine NATs
  nat_v_update <- nat_v

  who_any_v_time_since_dose <- (timestep - dose_time0[who_any_v]) * dt
  who_any_v_time_since_dose <- ceiling(who_any_v_time_since_dose)

  sub_mat <- cbind(who_any_v_time_since_dose, dose_num0[who_any_v])

  nat_v_update[who_any_v] <- nat_v[who_any_v] + (dr_vec_doses[sub_mat] * dt)

  # manual update: infection NATs
  nat_i_update <- nat_i

  who_any_i_time_since_dose <- (timestep - inf_time0[who_any_i]) * dt
  who_any_i_time_since_dose <- ceiling(who_any_i_time_since_dose)

  nat_i_update[who_any_i] <- nat_i[who_any_i] + (dr_vec_inf[who_any_i_time_since_dose] * dt)

  # safir update
  variables <- list(
    ab_titre = DoubleVariable$new(initial_values = nat_v),
    ab_titre_inf = DoubleVariable$new(initial_values = nat_i),
    dose_num = IntegerVariable$new(initial_values = dose_num0),
    dose_time = IntegerVariable$new(initial_values = dose_time0),
    inf_num = IntegerVariable$new(initial_values = inf_num0),
    inf_time = IntegerVariable$new(initial_values = inf_time0)
  )

  parameters <- list(
    dr_vec_doses = dr_vec_doses,
    dr_vec_inf = dr_vec_inf
  )

  proc <- independent_ab_titre_process(parameters = parameters, variables = variables, dt = dt)
  proc(timestep = timestep)
  for (v in variables) {
    v$.update()
  }

  nat_v_update_safir <- variables$ab_titre$get_values()
  nat_i_update_safir <- variables$ab_titre_inf$get_values()

  expect_equal(nat_v_update, nat_v_update_safir)
  expect_equal(nat_i_update, nat_i_update_safir)

})

