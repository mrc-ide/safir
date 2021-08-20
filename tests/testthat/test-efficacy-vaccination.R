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
  vaccinated$not(inplace = TRUE)

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

  # safir parameters
  parameters <- get_vaccine_ab_titre_parameters(vaccine = "Pfizer")
  parameters$N_phase <- 1

  n <- 1e2
  tmax <- 800

  parameters$population <- n

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

  dr_vec <- c(rep(dr_s, period_s),
              seq(dr_s, dr_l, length.out = time_to_decay),
              dr_l)

  mu_ab_d1 <- 0.14 # mean titre dose 1
  std10 <- 0.44 # Pooled standard deviation of antibody level on log10 data

  nt <- NULL
  t <- 0:tmax
  time_to_decay <- t_period_l - period_s # time in days to reach longest half-life

  # vector of decay rates over time: first we have a period of fast decay, then gradually shift to period of long decay
  dr_vec <- c(rep(dr_s, period_s),
              seq(dr_s, dr_l, length.out = time_to_decay),
              rep(dr_l, (length(t) - t_period_l)))

  # nab titre distribution for a given vaccine product: draw from a log-normal distribution

  z1 <- rnorm(n, log10(mu_ab_d1), std10)

  # initiate titre vector
  nt <- matrix(data = NaN,nrow = length(t),ncol = n)
  nt[1, ] <- log(10^z1)

  # decay antibodies over time on natural log scale
  for (i in (2:length(t))){
    nt[i, ] <- nt[i-1, ] + dr_vec[i]
  }

  # safir
  variables <- create_vaccine_variables(variables = list(),parameters = parameters)

  schedule_dose_vaccine(timestep = 0,variables = variables,target = Bitset$new(n)$insert(1:n),dose = 1,parameters = parameters)
  variables$ab_titre$queue_update(values = log(10^z1)) # make sure using the same set of RVs
  variables$ab_titre$.update()
  variables$dose_num$.update()
  variables$dose_time[[1]]$.update()

  ab_titre <- vaccine_ab_titre_process(parameters = parameters,variables = variables,events = list(),dt = 1)

  safir_out <- matrix(data = NaN,nrow = tmax + 1,ncol = n)
  safir_out[1, ] <- variables$ab_titre$get_values()
  for (t in 2:(tmax+1)) {
    ab_titre(timestep = t)
    variables$ab_titre$.update()
    safir_out[t, ] <-   variables$ab_titre$get_values()
  }

  # test Ab time series the same
  expect_equal(rowMeans(safir_out),rowMeans(nt))

  # test efficacy against infection is the same
  nt1 <- exp(nt[,1])
  ef_infection <- 1 / (1 + exp(-k * (log10(nt1) - log10(ab_50))))

  ef_infection_safir <- vaccine_efficacy_infection(ab_titre = safir_out[,1],parameters = parameters)
  ef_infection_safir <- 1 - ef_infection_safir

  expect_equal(ef_infection,ef_infection_safir)

  # test efficacy against severe disease is the same
  ef_severe_uncond <- 1 / (1 + exp(-k * (log10(nt1) - log10(ab_50_severe))))
  ef_severe <-  1 - ((1 - ef_severe_uncond)/(1 - ef_infection))

  ef_severe_safir <- vaccine_efficacy_severe(ab_titre = safir_out[,1],ef_infection = ef_infection_safir,parameters = parameters)
  ef_severe_safir <- 1 - ef_severe_safir

  expect_equal(ef_severe,ef_severe_safir)
})
