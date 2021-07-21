test_that("c++ infection process (multiple doses, no types) returns identical results as R version", {

  # pars
  iso3c <- "GBR"
  pop <- safir:::get_population(iso3c)
  contact_mat <- squire::get_mixing_matrix(iso3c = iso3c)
  pop$n <- as.integer(pop$n / 100)

  parameters <- safir::get_parameters(
    population = pop$n,
    contact_matrix_set = contact_mat,
    iso3c = iso3c,
    time_period = 365
  )
  parameters$beta_set <- parameters$beta_set*rexp(n = length(parameters$beta_set))

  # test states
  n <- 1e5
  dt <- 0.5
  valid_states <- c("S","IMild","ICase","IAsymp")
  state0 <- sample(x = valid_states,size = n,replace = T)
  age0 <- sample.int(n = 17,size = n,replace = T)
  ef_inf0 <- runif(n = n)

  # test R
  states <- individual::CategoricalVariable$new(categories = valid_states,initial_values = state0)
  discrete_age <- individual::IntegerVariable$new(initial_values = age0)
  ef_infection <- individual::DoubleVariable$new(initial_values = ef_inf0)
  exposure <- individual::TargetedEvent$new(population_size = n)

  set.seed(5436L)
  inf_proc_R <- infection_process_vaccine(parameters = parameters,variables = list(states=states,discrete_age=discrete_age,ef_infection=ef_infection),events = list(exposure=exposure),dt =dt)
  inf_proc_R(timestep = 100)

  sched_R <- exposure$get_scheduled()$to_vector()

  # test C++
  states1 <- individual::CategoricalVariable$new(categories = valid_states,initial_values = state0)
  discrete_age1 <- individual::IntegerVariable$new(initial_values = age0)
  ef_infection1 <- individual::DoubleVariable$new(initial_values = ef_inf0)
  exposure1 <- individual::TargetedEvent$new(population_size = n)

  set.seed(5436L)
  inf_proc_Cpp <- infection_process_vaccine_cpp(parameters = parameters,variables = list(states=states1,discrete_age=discrete_age1,ef_infection=ef_infection1),events = list(exposure=exposure1),dt =dt)
  individual:::execute_process(process = inf_proc_Cpp,timestep = 100)

  sched_Cpp <- exposure1$get_scheduled()$to_vector()

  expect_identical(sort(sched_R), sort(sched_Cpp))
})
