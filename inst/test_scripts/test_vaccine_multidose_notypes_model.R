# --------------------------------------------------------------------------------
#   test vaxx model, multi dose, no types
# --------------------------------------------------------------------------------

rm(list=ls());gc()
library(individual)

iso3c <- "GBR"
pop <- safir:::get_population(iso3c)
pop$n <- as.integer(pop$n / 1e3)
contact_mat <- squire::get_mixing_matrix(iso3c = iso3c)

tmax <- 10
dt <- 0.5
R0 <- 4

# vaccine dosing
vaccine_doses <- 2
dose_period <- c(NaN, 2)
vaccine_set <- c(0, rep(1e4, 9))

# vaccine strategy
vaccine_coverage_mat <- strategy_matrix(strategy = "Elderly",max_coverage = 0.2)
next_dose_priority <- matrix(data = 0, nrow = vaccine_doses - 1,ncol = ncol(vaccine_coverage_mat))
next_dose_priority[1, 15:17] <- 1 # prioritize 3 oldest age groups
storage.mode(next_dose_priority) <- "integer"

# base parameters
parameters <- safir::get_parameters(
  population = pop$n,
  contact_matrix_set = contact_mat,
  iso3c = iso3c,
  R0 = R0,
  time_period = tmax
)

# attach vaccine parameters
parameters$N_prioritisation_steps <- nrow(vaccine_coverage_mat)
parameters$vaccine_coverage_mat <- vaccine_coverage_mat
parameters$next_dose_priority <- next_dose_priority

parameters$vaccine_set <- vaccine_set
parameters$dose_period <- dose_period
parameters$N_phase <- vaccine_doses


# create variables
timesteps <- parameters$time_period/dt
variables <- create_variables(pop = pop, parameters = parameters)3
variables <- create_vaccine_variables(variables = variables,pop = pop$n,max_dose = vaccine_doses)

# create events
events <- create_events(parameters = parameters)
events <- create_events_vaccination(events = events,parameters = parameters)


safir::attach_event_listeners(variables = variables,events = events,parameters = parameters, dt = dt)
renderer <- individual::Render$new(timesteps)
processes <- list(
  safir::infection_process_cpp(parameters = parameters,variables = variables,events = events,dt = dt),
  individual::categorical_count_renderer_process(renderer, variables$state, categories = variables$states$get_categories())
)
safir::setup_events(parameters = parameters,events = events,variables = variables,dt = dt)

individual::simulation_loop(
  variables = variables,
  events = events,
  processes = processes,
  timesteps = timesteps
)
