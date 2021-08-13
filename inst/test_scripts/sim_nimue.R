rm(list=ls());gc();
# dev.off()
library(safir)
library(nimue)
library(individual)
library(data.table)
library(ggplot2)

iso3c <- "GBR"
pop <- safir:::get_population(iso3c)
pop$n <- as.integer(pop$n / 2e2)
contact_mat <- squire::get_mixing_matrix(iso3c = iso3c)

tmax <- 200
R0 <- 4

vaccine_coverage_mat <- strategy_matrix(strategy = "Elderly",max_coverage = 0.2)
tt_vaccine <- c(0, 10:100)
max_vaccine <- c(0, seq(1e3, 5e4, length.out = length(tt_vaccine)-1))

# nimue run
increasing <- run(
  time_period = tmax,
  population = pop$n,
  R0 = R0,
  contact_matrix_set = contact_mat,
  vaccine_coverage_mat = vaccine_coverage_mat,
  max_vaccine = max_vaccine,
  tt_vaccine = tt_vaccine
)

# safir run
dt <- 0.01

parameters <- get_parameters_nimue(
  population = pop$n,
  contact_mat = contact_mat,
  time_period = tmax,
  R0 = R0,
  max_vaccine = max_vaccine,
  tt_vaccine = tt_vaccine,
  vaccine_coverage_mat = vaccine_coverage_mat,
  dt = dt
)

timesteps <- parameters$time_period/dt

variables <- create_variables(pop = pop, parameters = parameters)
variables <- create_vaccine_variables_nimue(variables = variables,pop = pop)

events <- create_events(parameters = parameters)
events <- create_events_nimue(events = events,parameters = parameters)

attach_event_listeners(variables = variables,events = events,parameters = parameters, dt = dt)
attach_event_listeners_nimue(variables = variables,events = events,parameters = parameters,dt = dt)

renderer <- Render$new(parameters$time_period)
vaxx_renderer <- Render$new(parameters$time_period)

processes <- list(
  vaccination_process_nimue_cpp(parameters = parameters,variables = variables,events = events,dt = dt),
  infection_process_nimue_cpp(parameters = parameters,variables = variables,events = events,dt = dt),
  categorical_count_renderer_process_daily(renderer = renderer, variable = variables$state, categories = variables$states$get_categories(),dt = dt),
  integer_count_render_process_daily(renderer = vaxx_renderer,variable = variables$vaccine_states,margin = 1:4,dt = dt)
)

setup_events_nimue(parameters = parameters,events = events,variables = variables,dt = dt)

system.time(
  simulation_loop_nimue(
    variables = variables,
    events = events,
    processes = processes,
    timesteps = timesteps,
    progress = TRUE
  )
)

# extract data - state
nimue_compartments<- c("S","E","D","R","IMild","ICase","IRec","IICU","IHospital")

df <- renderer$to_dataframe()
safir_dt <- as.data.table(df)
safir_dt[, IMild_count := IMild_count + IAsymp_count]
safir_dt[, IAsymp_count := NULL]
safir_dt[, IICU_count := IMVNotGetDie_count + IMVNotGetLive_count + IMVGetLive_count + IMVGetDie_count]
safir_dt[, c("IMVNotGetDie_count","IMVNotGetLive_count","IMVGetLive_count","IMVGetDie_count") := NULL]
safir_dt[, IHospital := IOxGetDie_count + IOxNotGetDie_count + IOxNotGetLive_count + IOxGetLive_count]
safir_dt[, c("IOxGetDie_count", "IOxNotGetDie_count", "IOxNotGetLive_count", "IOxGetLive_count") := NULL]
safir_dt <- melt(safir_dt,id.vars = "timestep",variable.name = "compartment",value.name = "value")
safir_dt[, compartment := gsub("(^)(\\w*)(_count)", "\\2", compartment)]
setnames(safir_dt,old = "timestep",new = "t")
safir_dt <- safir_dt[compartment %in% nimue_compartments, ]
safir_dt[ ,"model" := "safir"]

nimue_dt <- as.data.table(format(increasing, compartments = nimue_compartments))
nimue_dt[ , replicate := NULL]
nimue_dt <- nimue_dt[compartment %in% nimue_compartments, ]
nimue_dt[ ,"model" := "nimue"]

combined_dt <- rbind(nimue_dt,safir_dt)

ggplot(data = combined_dt) +
  geom_line(aes(x=t,y=value,color=compartment,group=model,linetype=model)) +
  facet_wrap(.~compartment,scales="free_y")+
  guides(color = FALSE)+
  theme_bw()+
  theme(strip.text.x = element_text(size=10,face = "bold"))

# extract data - vaxx
safir_vax_dt <- as.data.table(vaxx_renderer$to_dataframe())
safir_vax_dt[ , "vaccinated" := X2_count + X3_count ]
safir_vax_dt[ , c("X2_count","X3_count") := NULL]
setnames(x = safir_vax_dt,old = c("X1_count","X4_count"),new = c("unvaccinated","priorvaccinated"))
safir_vax_dt <- melt(safir_vax_dt,id.vars = "timestep",variable.name = "compartment")
setnames(safir_vax_dt,old = "timestep",new = "t")
safir_vax_dt[, "model" := "safir"]

nimue_vax_dt <- as.data.table(format(increasing, compartments = NULL,summaries = c("unvaccinated","vaccinated","priorvaccinated")))
nimue_vax_dt[, replicate := NULL]
nimue_vax_dt[, "model" := "nimue"]

combined_vax_dt <- rbind(safir_vax_dt,nimue_vax_dt)

ggplot(data = combined_vax_dt) +
  geom_line(aes(x=t,y=value,color=compartment,linetype=model,group=model)) +
  facet_wrap(.~compartment,scales="free_y") +
  guides(color = FALSE)+
  theme_bw()+
  theme(strip.text.x = element_text(size=10,face = "bold"))
