test_that("c++ infection process (multiple doses, no types) returns identical results as R version", {

  library(nimue)

  # pars
  iso3c <- "GBR"
  pop <- safir::get_population(iso3c)
  contact_mat <- squire::get_mixing_matrix(iso3c = iso3c)
  pop$n <- as.integer(pop$n / 100)

  parameters <- safir::get_parameters(
    population = pop$n,
    contact_matrix_set = contact_mat,
    iso3c = iso3c,
    time_period = 365,
    dt = 1
  )
  parameters$beta_set <- parameters$beta_set*rexp(n = length(parameters$beta_set))

  # test states
  n <- 1e5
  dt <- 0.5
  valid_states <- c("S","IMild","ICase","IAsymp")
  state0 <- sample(x = valid_states,size = n,replace = T)
  age0 <- sample.int(n = 17,size = n,replace = T)
  vaccine_parameters <- get_vaccine_ab_titre_parameters(vaccine = "Pfizer")
  ab_titre0 <- log(10^rnorm(n = n, mean = log10(vaccine_parameters$mu_ab[1]),sd = vaccine_parameters$std10))
  parameters <- make_vaccine_parameters(safir_parameters = parameters,vaccine_ab_parameters = vaccine_parameters,vaccine_set = rep(100,365),dose_period = c(NaN, 10),strategy_matrix = nimue::strategy_matrix(strategy = "Elderly"),next_dose_priority_matrix = matrix(0,nrow = 1,ncol = 17))

  # test R
  states <- individual::CategoricalVariable$new(categories = valid_states,initial_values = state0)
  discrete_age <- individual::IntegerVariable$new(initial_values = age0)
  ab_titre <- individual::DoubleVariable$new(initial_values = ab_titre0)
  exposure <- individual::TargetedEvent$new(population_size = n)

  set.seed(5436L)
  inf_proc_R <- infection_process_vaccine(parameters = parameters,variables = list(states=states,discrete_age=discrete_age,ab_titre=ab_titre),events = list(exposure=exposure),dt =dt)
  inf_proc_R(timestep = 100)

  sched_R <- exposure$get_scheduled()$to_vector()

  # test C++
  states1 <- individual::CategoricalVariable$new(categories = valid_states,initial_values = state0)
  discrete_age1 <- individual::IntegerVariable$new(initial_values = age0)
  ab_titre1 <- individual::DoubleVariable$new(initial_values = ab_titre0)
  exposure1 <- individual::TargetedEvent$new(population_size = n)

  set.seed(5436L)
  inf_proc_Cpp <- infection_process_vaccine_cpp(parameters = parameters,variables = list(states=states1,discrete_age=discrete_age1,ab_titre=ab_titre1),events = list(exposure=exposure1),dt =dt)
  execute_process(process = inf_proc_Cpp,timestep = 100)

  sched_Cpp <- exposure1$get_scheduled()$to_vector()

  expect_identical(sort(sched_R), sort(sched_Cpp))
})
