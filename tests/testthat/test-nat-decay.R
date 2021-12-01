test_that("test NAT decay from vaccination works on dt = 1", {

  library(individual)

  n <- 10
  tmax <- 800

  # Alexandra's code
  hl_s <- 108 # Half life of antibody decay - short
  hl_l <- 3650 # Half life of antibody decay - long
  period_s <- 250
  t_period_l <- 365 # Time point at which to have switched to longest half-life
  time_to_decay <- t_period_l - period_s # time in days to reach longest half-life
  dr_s <- -log(2)/hl_s # Corresponding decay rate in days for half life above
  dr_l <- -log(2)/hl_l

  ab_50 <- 0.2 # titre relative to convalescent required to provide 50% protection from infection, on linear scale
  ab_50_severe <- 0.03
  k <- 2.94 # shape parameter of efficacy curve

  mu_ab_d1 <- 0.14 # mean titre dose 1
  std10 <- 0.44 # Pooled standard deviation of antibody level on log10 data

  nt <- NULL
  t <- 0:tmax
  time_to_decay <- t_period_l - period_s # time in days to reach longest half-life

  # vector of decay rates over time: first we have a period of fast decay, then gradually shift to period of long decay
  dr_vec <- c(rep(dr_s, period_s),
              seq(dr_s, dr_l, length.out = time_to_decay),
              rep(dr_l, (length(t) - t_period_l)))

  z1 <- rnorm(n, log10(mu_ab_d1), std10)

  # initiate titre vector
  nt <- matrix(data = NaN,nrow = length(t),ncol = n)
  nt[1, ] <- log(10^z1)

  # decay antibodies over time on natural log scale
  for (i in (2:length(t))){
    nt[i, ] <- nt[i-1, ] + dr_vec[i-1]
  }

  # safir
  parameters <- get_vaccine_ab_titre_parameters(vaccine = "Pfizer")
  parameters$population <- n
  parameters$N_phase <- 1
  parameters$dr_vec <- dr_vec

  variables <- create_vaccine_variables(variables = list(),parameters = parameters)

  dt <- 0.5
  timesteps <- tmax/dt

  schedule_dose_vaccine(timestep = 0,variables = variables,target = Bitset$new(n)$insert(1:n),dose = 1,parameters = parameters)
  variables$ab_titre$queue_update(values = log(10^z1)) # make sure using the same set of RVs
  variables$ab_titre$.update()
  variables$dose_num$.update()
  variables$dose_time$.update()

  ab_titre <- vaccine_ab_titre_process(parameters = parameters,variables = variables, dt = dt)

  safir_out <- matrix(data = NaN,nrow = timesteps + 1,ncol = n)
  safir_out[1, ] <- variables$ab_titre$get_values()

  for (t in 1:timesteps) {
    ab_titre(timestep = t)
    variables$ab_titre$.update()
    safir_out[t + 1L, ] <-   variables$ab_titre$get_values()
  }

  # test Ab time series the same
  expect_equal(safir_out[nrow(safir_out), ], nt[nrow(nt), ])

})


test_that("test NAT decay from vaccination works on dt = 0.5", {

  library(individual)

  n <- 10
  tmax <- 800

  # Alexandra's code
  hl_s <- 108 # Half life of antibody decay - short
  hl_l <- 3650 # Half life of antibody decay - long
  period_s <- 250
  t_period_l <- 365 # Time point at which to have switched to longest half-life
  time_to_decay <- t_period_l - period_s # time in days to reach longest half-life
  dr_s <- -log(2)/hl_s # Corresponding decay rate in days for half life above
  dr_l <- -log(2)/hl_l

  ab_50 <- 0.2 # titre relative to convalescent required to provide 50% protection from infection, on linear scale
  ab_50_severe <- 0.03
  k <- 2.94 # shape parameter of efficacy curve

  mu_ab_d1 <- 0.14 # mean titre dose 1
  std10 <- 0.44 # Pooled standard deviation of antibody level on log10 data

  nt <- NULL
  t <- 0:tmax
  time_to_decay <- t_period_l - period_s # time in days to reach longest half-life

  # vector of decay rates over time: first we have a period of fast decay, then gradually shift to period of long decay
  dr_vec <- c(rep(dr_s, period_s),
              seq(dr_s, dr_l, length.out = time_to_decay),
              rep(dr_l, (length(t) - t_period_l)))

  z1 <- rnorm(n, log10(mu_ab_d1), std10)

  # initiate titre vector
  nt <- matrix(data = NaN,nrow = length(t),ncol = n)
  nt[1, ] <- log(10^z1)

  # decay antibodies over time on natural log scale
  for (i in (2:length(t))){
    nt[i, ] <- nt[i-1, ] + dr_vec[i-1]
  }

  # safir
  parameters <- get_vaccine_ab_titre_parameters(vaccine = "Pfizer")
  parameters$population <- n
  parameters$N_phase <- 1
  parameters$dr_vec <- dr_vec

  variables <- create_vaccine_variables(variables = list(),parameters = parameters)

  dt <- 0.5
  timesteps <- tmax/dt

  schedule_dose_vaccine(timestep = 0,variables = variables,target = Bitset$new(n)$insert(1:n),dose = 1,parameters = parameters)
  variables$ab_titre$queue_update(values = log(10^z1)) # make sure using the same set of RVs
  variables$ab_titre$.update()
  variables$dose_num$.update()
  variables$dose_time$.update()

  ab_titre <- vaccine_ab_titre_process(parameters = parameters,variables = variables, dt = dt)

  safir_out <- matrix(data = NaN,nrow = timesteps + 1,ncol = n)
  safir_out[1, ] <- variables$ab_titre$get_values()

  for (t in 1:timesteps) {
    ab_titre(timestep = t)
    variables$ab_titre$.update()
    safir_out[t + 1L, ] <-   variables$ab_titre$get_values()
  }

  # test Ab time series the same
  expect_equal(safir_out[nrow(safir_out), ], nt[nrow(nt), ])

})


test_that("test NAT decay from infection works on dt = 0.5", {

  library(individual)

  n <- 10
  tmax <- 800

  # Alexandra's code
  hl_s <- 108 # Half life of antibody decay - short
  hl_l <- 3650 # Half life of antibody decay - long
  period_s <- 250
  t_period_l <- 365 # Time point at which to have switched to longest half-life
  time_to_decay <- t_period_l - period_s # time in days to reach longest half-life
  dr_s <- -log(2)/hl_s # Corresponding decay rate in days for half life above
  dr_l <- -log(2)/hl_l

  ab_50 <- 0.2 # titre relative to convalescent required to provide 50% protection from infection, on linear scale
  ab_50_severe <- 0.03
  k <- 2.94 # shape parameter of efficacy curve

  mu_ab_d1 <- 0.14 # mean titre dose 1
  std10 <- 0.44 # Pooled standard deviation of antibody level on log10 data

  nt <- NULL
  t <- 0:tmax
  time_to_decay <- t_period_l - period_s # time in days to reach longest half-life

  # vector of decay rates over time: first we have a period of fast decay, then gradually shift to period of long decay
  dr_vec <- c(rep(dr_s, period_s),
              seq(dr_s, dr_l, length.out = time_to_decay),
              rep(dr_l, (length(t) - t_period_l)))

  z1 <- rnorm(n, log10(mu_ab_d1), std10)

  # initiate titre vector
  nt <- matrix(data = NaN,nrow = length(t),ncol = n)
  nt[1, ] <- log(10^z1)

  # decay antibodies over time on natural log scale
  for (i in (2:length(t))){
    nt[i, ] <- nt[i-1, ] + dr_vec[i-1]
  }

  # safir
  parameters <- get_vaccine_ab_titre_parameters(vaccine = "Pfizer")
  parameters$population <- n
  parameters$N_phase <- 1
  parameters$dr_vec <- dr_vec

  variables <- create_vaccine_variables(variables = list(),parameters = parameters)
  variables <- create_natural_immunity_variables(variables = variables, parameters = parameters)

  dt <- 0.5
  timesteps <- tmax/dt

  target <- Bitset$new(n)$insert(1:n)
  variables$inf_num$queue_update(values = 1, index = target)
  variables$inf_time$queue_update(values = 0, index = target)
  variables$ab_titre$queue_update(values = log(10^z1)) # make sure using the same set of RVs
  variables$ab_titre$.update()
  variables$inf_num$.update()
  variables$inf_time$.update()

  ab_titre <- natural_immunity_ab_titre_process(parameters = parameters,variables = variables, dt = dt)

  safir_out <- matrix(data = NaN,nrow = timesteps + 1,ncol = n)
  safir_out[1, ] <- variables$ab_titre$get_values()

  for (t in 1:timesteps) {
    ab_titre(timestep = t)
    variables$ab_titre$.update()
    safir_out[t + 1L, ] <-   variables$ab_titre$get_values()
  }

  # test Ab time series the same
  expect_equal(safir_out[nrow(safir_out), ], nt[nrow(nt), ])

})

