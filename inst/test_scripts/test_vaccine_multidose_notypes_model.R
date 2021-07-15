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
variables <- create_variables(pop = pop, parameters = parameters)
variables <- create_vaccine_variables(variables = variables,pop = pop$n,max_dose = vaccine_doses)

# create events
events <- create_events(parameters = parameters)
events <- create_events_vaccination(events = events,parameters = parameters)
attach_event_listeners(variables = variables,events = events,parameters = parameters, dt = dt)
attach_event_listeners_vaccination(variables = variables,events = events,parameters = parameters,dt = dt)

# make renderer
renderer <- Render$new(timesteps)

# processes
processes <- list(
  infection_process_vaccine(parameters = parameters,variables = variables,events = events,dt = dt),
  vaccine_ab_titre_process(parameters = parameters,variables = variables,events = events,dt = dt),
  vaccination_process(parameters = parameters,variables = variables,events = events,dt = dt),
  categorical_count_renderer_process(renderer, variables$state, categories = variables$states$get_categories())
)

setup_events_vaccine(parameters = parameters,events = events,variables = variables,dt = dt)

# if starting with people already vaccinated, use this function
# initialize_vaccine_variables(variables = ,dose_time_init = ,dose_num_init = )

# debug all processes
invisible(sapply(processes,debug))
debug(events$exposure$.listeners[[2]])
invisible(sapply(events$scheduled_dose,function(e){debug(e$.listeners[[1]])}))

simulation_loop_vaccine(
  variables = variables,
  events = events,
  processes = processes,
  timesteps = timesteps
)
