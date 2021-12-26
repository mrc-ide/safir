rm(list=ls());gc()
library(safir)
library(individual)
library(data.table)
library(ggplot2)
library(parallel)

iso3c <- "GBR"
pop <- safir:::get_population(iso3c)
contact_mat <- squire::get_mixing_matrix(iso3c = iso3c)

nrep <- 1
# Scale it for speed
pop$n <- as.integer(pop$n / 100)

# Create our simulation parameters
R0 <- 2
time_period <- 200
dt <- 0.1

parameters <- safir::get_parameters(
  population = pop$n,
  contact_matrix_set = contact_mat,
  iso3c = iso3c,
  R0 = R0,
  time_period = time_period,
  dt = dt
)


out <- squire::run_explicit_SEEIR_model(
  population = pop$n,
  country = "United Kingdom",
  contact_matrix_set = contact_mat,
  time_period = 200,
  replicates = nrep,
  day_return = TRUE,
  R0 = R0,
  dt = dt
)

timesteps <- parameters$time_period/dt
variables <- safir::create_variables(pop = pop, parameters = parameters)
events <- safir::create_events(parameters = parameters)
safir::attach_event_listeners(variables = variables,events = events,parameters = parameters, dt = dt)
renderer <- individual::Render$new(timesteps)
# hosp_render <- create_hosp_renderers(parameters)
# attach_hosp_listeners(renderers = hosp_render, events = events)
processes <- list(
  safir::infection_process_cpp(parameters = parameters,variables = variables,events = events,dt = dt),
  individual::categorical_count_renderer_process(renderer, variables$state, categories = variables$states$get_categories())
)
safir::setup_events(parameters = parameters,events = events,variables = variables,dt = dt)

system.time(individual::simulation_loop(
  variables = variables,
  events = events,
  processes = processes,
  timesteps = timesteps
))
df <- renderer$to_dataframe()

state_t0 <- get_state_vector(psq = parameters)

saf_dt <- as.data.table(df)
saf_dt <- rbind(data.table(t(state_t0)), saf_dt)
saf_dt[, IMild_count := IMild_count + IAsymp_count]
saf_dt[, IAsymp_count := NULL]
saf_dt <- melt(saf_dt,id.vars = c("timestep"),variable.name = "name")
saf_dt[, model := "safir"]
saf_dt[, name := gsub("(^)(\\w*)(_count)", "\\2", name)]
setnames(x = saf_dt,old = c("timestep","name","value"),new = c("t","compartment","y"))
saf_dt[, t := t * dt]

sq_dt <- as.data.table(squire::format_output(out, unique(saf_dt$compartment)))
sq_dt[, model := "squire"]
sq_dt[, replicate := NULL]

ggplot(data = rbind(saf_dt,sq_dt), aes(t,y,color = model)) +
  geom_line() +
  geom_line() +
  facet_wrap(~compartment, scales = "free")

