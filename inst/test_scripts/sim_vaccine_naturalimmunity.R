# --------------------------------------------------------------------------------
#   model with ab dynamics and natural immunity
# --------------------------------------------------------------------------------

rm(list=ls());gc()
library(individual)
library(data.table)
library(ggplot2)
library(safir)
library(nimue)

iso3c <- "GBR"
pop <- safir:::get_population(iso3c)
pop$n <- as.integer(pop$n / 1e3)
contact_mat <- squire::get_mixing_matrix(iso3c = iso3c)

tmax <- 100
dt <- 0.5
R0 <- 4

# vaccine dosing
vaccine_doses <- 2
dose_period <- c(NaN, 28)
vaccine_set <- vaccine_set <- c(0, seq(from = 1e3, to = 1e4, length.out = tmax-1))
vaccine_set <- floor(vaccine_set)

# vaccine strategy
vaccine_coverage_mat <- strategy_matrix(strategy = "Elderly",max_coverage = 0.5)
next_dose_priority <- matrix(data = 0, nrow = vaccine_doses - 1,ncol = ncol(vaccine_coverage_mat))
next_dose_priority[1, 15:17] <- 1 # prioritize 3 oldest age groups for next dose

# base parameters
parameters <- safir::get_parameters(
  population = pop$n,
  contact_matrix_set = contact_mat,
  iso3c = iso3c,
  R0 = R0,
  time_period = tmax,
  dt = dt,
  seeding_cases = 50
)

# vaccine parameters
ab_parameters <- get_vaccine_ab_titre_parameters(vaccine = "Pfizer", max_dose = vaccine_doses,correlated = TRUE)

# combine parameters and verify
parameters <- make_vaccine_parameters(
  safir_parameters = parameters,
  vaccine_ab_parameters = ab_parameters,
  vaccine_set = vaccine_set,
  dose_period = dose_period,
  strategy_matrix = vaccine_coverage_mat,
  next_dose_priority_matrix = next_dose_priority
)

# ab boost for each infection...assumes someone will not get infected more than 51 times.
parameters$mu_ab_infection <- rep(ab_parameters$mu_ab, times = c(1, 50))

# create variables
timesteps <- parameters$time_period/dt
variables <- create_variables(pop = pop, parameters = parameters)
variables <- create_vaccine_variables(variables = variables,parameters = parameters)
variables <- create_natural_immunity_variables(variables = variables, parameters = parameters)

# create events
events <- create_events(parameters = parameters)
events <- create_events_vaccination(events = events,parameters = parameters)
attach_event_listeners(variables = variables,events = events,parameters = parameters, dt = dt)
attach_event_listeners_vaccination(variables = variables,events = events,parameters = parameters,dt = dt)
attach_event_listeners_natural_immunity(variables = variables, events = events, parameters = parameters, dt = dt)

# make renderers
renderer <- Render$new(parameters$time_period)
ab_renderer <- matrix(data = NaN,nrow = parameters$time_period,ncol = sum(parameters$population))
dose_renderer <- Render$new(parameters$time_period)

double_count_render_process_daily <- function(variable, dt) {
  stopifnot(inherits(variable, "DoubleVariable"))
  function(t) {
    if ((t * dt) %% 1 == 0) {
      ab_renderer[as.integer(t * dt), ] <<- variable$get_values()
    }
  }
}

# processes
processes <- list(
  natural_immunity_ab_titre_process(parameters = parameters,variables = variables,events = events,dt = dt),
  vaccination_process(parameters = parameters,variables = variables,events = events,dt = dt),
  infection_process_vaccine_cpp(parameters = parameters,variables = variables,events = events,dt = dt),
  categorical_count_renderer_process_daily(renderer = renderer,variable = variables$states,categories = variables$states$get_categories(),dt = dt),
  double_count_render_process_daily(variable = variables$ab_titre,dt = dt),
  integer_count_render_process_daily(renderer = dose_renderer,variable = variables$dose_num,margin = 0:vaccine_doses,dt = dt)
)

setup_events_vaccine(parameters = parameters,events = events,variables = variables,dt = dt)

# stuff we'd like to check
debug(processes[[1]])
debugonce(events$recovery$.listeners[[2]])
debug(events$recovery$.listeners[[3]])

system.time(simulation_loop_safir(
  variables = variables,
  events = events,
  processes = processes,
  timesteps = timesteps,
  variables_dont_update = c("discrete_age", "phase"),
  progress = TRUE
))

