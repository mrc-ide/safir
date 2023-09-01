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

test_that("efficacy functions work in R and C++", {

  n <- 50
  vaccine_parameters <- get_vaccine_ab_titre_parameters(vaccine = "Pfizer")
  ab_titre <- log(10^rnorm(n = n, mean = log10(vaccine_parameters$mu_ab[1]),sd = vaccine_parameters$std10))
  ab_titre[sample.int(n = n,size = 20,replace = FALSE)] <- -Inf
  nat <- exp(ab_titre)

  ef_transmission_R <- vaccine_efficacy_infection(nat = nat, parameters = vaccine_parameters, day = 1)
  ef_transmission_cpp <- vaccine_efficacy_infection_cpp(nat = nat, parameters = vaccine_parameters, day = 1)

  expect_equal(ef_transmission_R, ef_transmission_cpp)

})


test_that("natural_immunity_ab_titre_process works with mixed inf/vaccinated", {

  N_phase <- 3

  # base parameters
  parameters <- safir::get_parameters(
    population = rep(1e2, 17),
    iso3c = "GBR",
    R0 = 4,
    time_period = 10,
    dt = 1
  )

  # vaccine parameters
  mu_ab_list <- data.frame(name = c("Pfizer", "AstraZeneca", "Sinovac", "Moderna"),
                           mu_ab_d1 = c(13/94, 1/59, 28/164, ((185+273)/2)/321),
                           mu_ab_d2 = c(223/94, 32/59, 28/164, 654/158),
                           mu_ab_d3 = c(223/94, 32/59, 28/164, 654/158))

  ab_parameters <- get_vaccine_ab_titre_parameters(vaccine = "Pfizer", max_dose = N_phase,correlated = FALSE, mu_ab_list = mu_ab_list)

  n <- 1e3
  dt <- 0.5
  dose_num_vec <- sample(x = 0:3, size = n, replace = TRUE)
  inf_num_vec <- sample(x = 0:2,size = n,replace = TRUE)
  who <- vector("list",3)
  when <- vector("list",3)

  dose_num <- IntegerVariable$new(initial_values = dose_num_vec)

  dose_time_vec <- rep(0, n)
  dose_windows <- list(c(1,20), c(21, 30), c(31, 40))
  for (i in 1:3) {
    dosed <- which(dose_num_vec == i)
    if (length(dosed) > 0) {
      dose_time_vec[dosed] <- sample(x = dose_windows[[i]][1]:dose_windows[[i]][2], size = length(dosed), replace = TRUE)
    }
  }
  dose_time <- IntegerVariable$new(initial_values = dose_time_vec)

  inf_num <- IntegerVariable$new(initial_values = inf_num_vec)
  inf_time <- IntegerVariable$new(initial_values = rep(-1, n))
  inf_time$queue_update(values = sample(x = 1:10, size = sum(inf_num_vec == 1), replace = TRUE), index = which(inf_num_vec == 1))
  inf_time$queue_update(values = sample(x = 10:20, size = sum(inf_num_vec == 2), replace = TRUE), index = which(inf_num_vec == 2))
  inf_time$.update()

  vaccinated <- dose_num$get_index_of(set = 0)$not(inplace = TRUE)
  infected <- inf_num$get_index_of(set = 0)$not(inplace = TRUE)

  vaccinated_or_infected <- vaccinated$or(infected)

  ab_0 <- rep(-Inf, n)
  ab_0[vaccinated_or_infected$to_vector()] <- 1
  ab_titre <- DoubleVariable$new(initial_values = ab_0)

  # safir
  vars <- list(ab_titre = ab_titre, dose_num = dose_num, dose_time = dose_time, inf_num = inf_num, inf_time = inf_time)
  proc <- natural_immunity_ab_titre_process(parameters = c(parameters, ab_parameters), variables = vars, dt = dt)
  proc(timestep = 41)

  vars$ab_titre$.update()
  ab_after_update <- vars$ab_titre$get_values()

  expect_true(all(ab_0[vaccinated_or_infected$to_vector()] != ab_after_update[vaccinated_or_infected$to_vector()]))
  expect_equal(sum(is.finite(ab_after_update)), vaccinated_or_infected$size())

})
