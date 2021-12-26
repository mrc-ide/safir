test_that("cross tab for compartments/age works", {

  SIR <- c("S", "I", "R")
  parameters <- list(N_age = 6)

  ages <- sample.int(n = parameters$N_age, size = 100, replace = TRUE)
  compartments <- sample(x = SIR, size = 100, replace = TRUE)
  sir_to_int <- c("S" = 1, "I" = 2, "R" = 3)

  comp1 <- matrix(data = 0, nrow = parameters$N_age, ncol = 3)
  for (i in 1:100) {
    comp1[ages[i], sir_to_int[[compartments[i]]]] <- comp1[ages[i], sir_to_int[[compartments[i]]]] + 1
  }

  expect_equal(colSums(comp1), as.vector(table(compartments)[SIR]))
  expect_equal(rowSums(comp1), as.vector(table(ages)[as.character(1:parameters$N_age)]))

  comp_variable <- CategoricalVariable$new(categories = SIR, initial_values = compartments)
  age_variable <- IntegerVariable$new(initial_values = ages)

  comp_render <- Render$new(timesteps = 1)

  render_proc <- compartments_age_render_process_daily(renderer = comp_render, age = age_variable, compartments = comp_variable, parameters = parameters, dt = 1)
  render_proc(t = 0.25) # no effect unless it hits one day
  render_out <- comp_render$to_dataframe()

  expect_equal(ncol(render_out), 1)

  render_proc(t = 1)
  render_out <- comp_render$to_dataframe()
  render_out <- render_out[, -1]

  comp2 <- matrix(data = 0, nrow = parameters$N_age, ncol = 3)

  for (i in 1:ncol(render_out)) {
    c_i <- strsplit(x = colnames(render_out)[i], split = "_")[[1]][2]
    a_i <- as.integer(strsplit(x = colnames(render_out)[i], split = "_")[[1]][4])
    comp2[a_i , sir_to_int[[c_i]]] <- render_out[1, i]
  }

  expect_equal(comp1, comp2)
})


test_that("hospitalization/ICU rendering creation works", {

  parameters <- get_parameters(iso3c = "GBR", dt = 0.5, time_period = 5)
  hosp_render <- create_hosp_renderers(parameters = parameters)

  expect_true(length(hosp_render) == 8L)
  expect_true(all(names(hosp_render) %in% c("ICU_get_live", "ICU_get_die", "hosp_get_live", "hosp_get_die", "ICU_not_get_live", "ICU_not_get_die", "hosp_not_get_live", "hosp_not_get_die")))
  expect_true(all(vapply(X = hosp_render, FUN = function(x) {inherits(x, "Render")}, FUN.VALUE = logical(1))))
})


test_that("incidence + hosp/ICU output is working correctly", {

  R0 <- 15
  time_period <- 2
  dt <- 1

  pop <- get_population("GBR")
  pop$n <- rep(1e3, 17) + rep(50, 17)

  parameters <- safir::get_parameters(
    population = pop$n,
    seeding_cases = 0,
    iso3c = "GBR",
    R0 = R0,
    time_period = time_period,
    dt = dt
  )
  parameters$S_0 <- parameters$S_0 - rep(50, 17)
  parameters$ICase1_0 <- rep(50, 17)

  timesteps <- parameters$time_period/dt
  variables <- create_variables(pop = pop, parameters = parameters)
  events <- create_events(parameters = parameters)
  attach_event_listeners(variables = variables,events = events,parameters = parameters, dt = dt)
  renderer <- individual::Render$new(timesteps)
  hosp_render <- create_hosp_renderers(parameters = parameters)
  attach_hosp_listeners(renderers = hosp_render, events = events)
  processes <- list(
    infection_process_cpp(parameters = parameters,variables = variables,events = events,dt = dt),
    individual::categorical_count_renderer_process(renderer, variables$state, categories = variables$states$get_categories())
  )
  setup_events(parameters = parameters,events = events,variables = variables,dt = dt)

  # incidence render
  incidence_renderer <- individual::Render$new(timesteps)
  attach_tracking_listener_incidence(events = events, renderer = incidence_renderer)

  # age-stratified incidence render
  age_incidence_renderer <- individual::Render$new(timesteps)
  attach_tracking_listener_age_incidence(events = events, renderer = age_incidence_renderer, age = variables$discrete_age, parameters = parameters)

  individual::simulation_loop(
    variables = variables,
    events = events,
    processes = processes,
    timesteps = timesteps
  )

  out_compartment <- renderer$to_dataframe()
  out_inc <- incidence_renderer$to_dataframe()
  out_age_inc <- age_incidence_renderer$to_dataframe()

  expect_true(all.equal(out_compartment[2, "E_count"], out_inc[1, "incidence"], sum(out_age_inc[1, -1])))

  out_hosp <- process_hosp_renderers(renderers = hosp_render, parameters = parameters)
  expect_true(nrow(out_hosp) == 2L)
})

