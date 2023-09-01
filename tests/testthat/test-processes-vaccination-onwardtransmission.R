test_that("c++ infection process (multiple doses, no types) testing with NAT affecting onward transmission", {

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
  ab_titre0 <- rep(-Inf, n)
  parameters <- make_vaccine_parameters(safir_parameters = parameters,vaccine_ab_parameters = vaccine_parameters,vaccine_set = rep(100,365),dose_period = c(NaN, 10),strategy_matrix = nimue::strategy_matrix(strategy = "Elderly"),next_dose_priority_matrix = matrix(0,nrow = 1,ncol = 17))

  # standard (no NAT effect)
  states <- individual::CategoricalVariable$new(categories = valid_states,initial_values = state0)
  discrete_age <- individual::IntegerVariable$new(initial_values = age0)
  ab_titre <- individual::DoubleVariable$new(initial_values = ab_titre0)
  ab_titre_inf <- individual::DoubleVariable$new(initial_values = ab_titre0)
  exposure <- individual::TargetedEvent$new(population_size = n)

  vars <- list(states=states,discrete_age=discrete_age,ab_titre=ab_titre,ab_titre_inf=ab_titre_inf)

  set.seed(5436L)
  inf_proc_no_nat <- infection_process_vaccine(parameters = parameters,variables = vars,events = list(exposure=exposure),dt = dt)
  inf_proc_no_nat(timestep = 100)

  no_nat <- exposure$get_scheduled()$to_vector()

  # NAT effect (with no titre)
  parameters$nt_efficacy_transmission <- TRUE
  exposure <- individual::TargetedEvent$new(population_size = n)

  set.seed(5436L)
  inf_proc_nat <- infection_process_vaccine(parameters = parameters,variables = vars,events = list(exposure=exposure),dt = dt)
  inf_proc_nat(timestep = 100)

  nat_no_titre <- exposure$get_scheduled()$to_vector()

  expect_equal(no_nat, nat_no_titre)

  # make sure significantly less infections occur when using NAT effect with a large ab titre
  zdose <- log(10^rnorm(n = n, mean = log10(3), sd = parameters$std10))
  ab_titre <- individual::DoubleVariable$new(initial_values = zdose)
  ab_titre_inf <- individual::DoubleVariable$new(initial_values = ab_titre0)
  exposure <- individual::TargetedEvent$new(population_size = n)

  vars <- list(states=states,discrete_age=discrete_age,ab_titre=ab_titre,ab_titre_inf=ab_titre_inf)

  set.seed(5436L)
  inf_proc_nat <- infection_process_vaccine(parameters = parameters,variables = vars,events = list(exposure=exposure),dt = dt)
  inf_proc_nat(timestep = 100)

  nat_high_titre <- exposure$get_scheduled()$to_vector()

  expect_true(length(nat_high_titre) < length(nat_no_titre))

})


test_that("R/C++ infection process consistent with NAT onward infectiousness and random draws of ab titre", {

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
  parameters <- make_vaccine_parameters(safir_parameters = parameters,vaccine_ab_parameters = vaccine_parameters,vaccine_set = rep(100,365),dose_period = c(NaN, 10),strategy_matrix = nimue::strategy_matrix(strategy = "Elderly"),next_dose_priority_matrix = matrix(0,nrow = 1,ncol = 17))
  parameters$nt_efficacy_transmission <- TRUE

  ab_titre0 <- log(10^rnorm(n = n, mean = log10(1), sd = parameters$std10))

  # R
  states <- individual::CategoricalVariable$new(categories = valid_states,initial_values = state0)
  discrete_age <- individual::IntegerVariable$new(initial_values = age0)
  ab_titre <- individual::DoubleVariable$new(initial_values = ab_titre0)
  ab_titre_inf <- individual::DoubleVariable$new(initial_values = ab_titre0)
  exposure <- individual::TargetedEvent$new(population_size = n)

  vars <- list(states=states,discrete_age=discrete_age,ab_titre=ab_titre,ab_titre_inf=ab_titre_inf)

  set.seed(1967391L)
  inf_proc_nat <- infection_process_vaccine(parameters = parameters,variables = vars,events = list(exposure=exposure),dt = dt)
  inf_proc_nat(timestep = 100)

  inf_R <- exposure$get_scheduled()$to_vector()

  # C++
  exposure <- individual::TargetedEvent$new(population_size = n)
  set.seed(1967391L)
  inf_proc_nat <- infection_process_vaccine_cpp(parameters = parameters,variables = vars,events = list(exposure=exposure),dt = dt)
  execute_process(process = inf_proc_nat,timestep = 100)

  inf_cpp <- exposure$get_scheduled()$to_vector()

  expect_equal(inf_R, inf_cpp)

})

test_that("R/C++ infection process consistent with NAT onward infectiousness and zero (-Inf) ab titre", {

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
  parameters <- make_vaccine_parameters(safir_parameters = parameters,vaccine_ab_parameters = vaccine_parameters,vaccine_set = rep(100,365),dose_period = c(NaN, 10),strategy_matrix = nimue::strategy_matrix(strategy = "Elderly"),next_dose_priority_matrix = matrix(0,nrow = 1,ncol = 17))
  parameters$nt_efficacy_transmission <- TRUE

  ab_titre0 <- rep(-Inf, n)

  # R
  states <- individual::CategoricalVariable$new(categories = valid_states,initial_values = state0)
  discrete_age <- individual::IntegerVariable$new(initial_values = age0)
  ab_titre <- individual::DoubleVariable$new(initial_values = ab_titre0)
  ab_titre_inf <- individual::DoubleVariable$new(initial_values = ab_titre0)
  exposure <- individual::TargetedEvent$new(population_size = n)

  vars <- list(states=states,discrete_age=discrete_age,ab_titre=ab_titre,ab_titre_inf=ab_titre_inf)

  set.seed(1967391L)
  inf_proc_nat <- infection_process_vaccine(parameters = parameters,variables = vars,events = list(exposure=exposure),dt = dt)
  inf_proc_nat(timestep = 100)

  inf_R <- exposure$get_scheduled()$to_vector()

  # C++
  exposure <- individual::TargetedEvent$new(population_size = n)
  set.seed(1967391L)
  inf_proc_nat <- infection_process_vaccine_cpp(parameters = parameters,variables = vars,events = list(exposure=exposure),dt = dt)
  execute_process(process = inf_proc_nat,timestep = 100)

  inf_cpp <- exposure$get_scheduled()$to_vector()

  expect_equal(inf_R, inf_cpp)

})

