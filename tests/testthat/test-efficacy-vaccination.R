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


test_that("vaccine_ab_titre_process works for everyone on dose 1", {

  hl_s <- 108 # Half life of antibody decay - short
  hl_l <- 3650 # Half life of antibody decay - long
  period_s <- 250
  t_period_l <- 365 # Time point at which to have switched to longest half-life
  time_to_decay <- t_period_l - period_s # time in days to reach longest half-life
  dr_s <- -log(2)/hl_s # Corresponding decay rate in days for half life above
  dr_l <- -log(2)/hl_l

  dr_vec <- c(rep(dr_s, period_s),
              seq(dr_s, dr_l, length.out = time_to_decay),
              dr_l)

  mu_ab_d1 <- 0.14 # mean titre dose 1
  std10 <- 0.44 # Pooled standard deviation of antibody level on log10 data

  parameters <- list(
    dr_vec = dr_vec,
    mu_ab_d1 = mu_ab_d1,
    std10 = std10
  )

  variables <- create_vaccine_variables(variables = list(),pop = 100,max_dose = 1)

})
