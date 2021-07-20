rm(list=ls());gc()
library(safir)
library(data.table)
library(ggplot2)
library(parallel)

iso3c <- "GBR"
pop <- safir:::get_population(iso3c)
contact_mat <- squire::get_mixing_matrix(iso3c = iso3c)

# use as many as you want normally.
options("mc.cores" = 20)

nrep <- 20
# Scale it for speed
pop$n <- as.integer(pop$n / 100)

# Create our simulation parameters
R0 <- 2
time_period <- 200

parameters <- safir::get_parameters(
  population = pop$n,
  contact_matrix_set = contact_mat,
  iso3c = iso3c,
  R0 = R0,
  time_period = time_period
)

dt <- 0.1


out <- squire::run_explicit_SEEIR_model(
  population = pop$n,
  country = "United Kingdom",
  contact_matrix_set = contact_mat,
  time_period = 200,
  replicates = nrep,
  day_return = TRUE,
  R0 = 2,
  dt = dt
)


system.time(
  saf_reps <- mclapply(X = 1:nrep,FUN = function(x){

    timesteps <- parameters$time_period/dt
    variables <- safir::create_variables(pop = pop, parameters = parameters)
    events <- safir::create_events(parameters = parameters)
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
    df <- renderer$to_dataframe()
    df$repetition <- x
    return(df)
  })
)


saf_reps <- do.call(rbind,saf_reps)

saf_dt <- as.data.table(saf_reps)
saf_dt[, IMild_count := IMild_count + IAsymp_count]
saf_dt[, IAsymp_count := NULL]
saf_dt <- melt(saf_dt,id.vars = c("timestep","repetition"),variable.name = "name")
saf_dt[, model := "safir"]
saf_dt[, name := gsub("(^)(\\w*)(_count)", "\\2", name)]
setnames(x = saf_dt,old = c("timestep","name","value"),new = c("t","compartment","y"))
saf_dt[, t := t * dt]
saf_dt <- saf_dt[, .(ymin = quantile(y,0.025), ymax = quantile(y,0.975), y = median(y)), by = .(t,compartment,model)]


sq_dt <- as.data.table(squire::format_output(out, unique(saf_dt$compartment)))
sq_dt[, model := "squire"]
sq_dt <- sq_dt[, .(ymin = quantile(y,0.025), ymax = quantile(y,0.975), y = median(y)), by = .(t,compartment,model)]


ggplot(data = rbind(saf_dt,sq_dt), aes(t,y,color = model)) +
  geom_line() +
  geom_ribbon(ggplot2::aes(ymin = ymin, ymax = ymax, fill = model), alpha = 0.2) +
  geom_line() +
  facet_wrap(~compartment, scales = "free")

