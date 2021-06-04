rm(list=ls());gc();dev.off()
library(safir)
library(nimue)
library(individual)
library(data.table)
library(ggplot2)

iso3c <- "GBR"
pop <- safir:::get_population(iso3c)
pop$n <- as.integer(pop$n / 2e2)
contact_mat <- squire::get_mixing_matrix(iso3c = iso3c)

tmax <- 365
R0 <- 4

# nimue run
increasing <- run(
  time_period = tmax,
  population = pop$n,
  R0 = R0,
  contact_matrix_set = contact_mat,
  max_vaccine = c(0, seq(10, 1e4, length.out = 40)),
  tt_vaccine = c(0, seq(10, 50, length.out = 40)),
  vaccine_efficacy_disease = rep(0, 17),
  vaccine_efficacy_infection = rep(0.9, 17)
)


# safir runs
dt <- 0.1
nrep <- 40
options("mc.cores" = 20)

parameters <- safir::get_parameters(
  population = pop$n,
  contact_matrix_set = contact_mat,
  iso3c = iso3c,
  R0 = R0,
  time_period = tmax
)

parameters <- append_vaccine_nimue(
  parameters = parameters,
  population = pop$n,
  contact_mat = contact_mat,
  max_vaccine = c(0, seq(10, 1e4, length.out = 40)),
  tt_vaccine = c(0, seq(10, 50, length.out = 40)),
  vaccine_efficacy_disease = rep(0, 17),
  vaccine_efficacy_infection = rep(0.9, 17)
)

timesteps <- parameters$time_period/dt

system.time(
  safir_reps <- mclapply(X = 1:nrep,FUN = function(x){

    variables <- create_variables(pop = pop, parameters = parameters)
    variables <- create_vaccine_variables_nimue(variables = variables,pop = pop)

    events <- create_events(parameters = parameters)
    events <- create_events_nimue(events = events,parameters = parameters)

    attach_event_listeners(variables = variables,events = events,parameters = parameters, dt = dt)
    attach_event_listeners_nimue(variables = variables,events = events,parameters = parameters,dt = dt)

    # this is bad and i should feel bad
    events$exposure$.listeners[[2]] <- NULL
    events$exposure$add_listener(
      safir:::create_exposure_scheduler_listener_nimue(
        events = events,
        variables = variables,
        parameters = parameters,
        dt = dt
      )
    )

    renderer <- Render$new(timesteps)
    processes <- list(
      vaccination_process_nimue(parameters = parameters,variables = variables,events = events,dt = dt),
      infection_process_nimue(parameters = parameters,variables = variables,events = events,dt = dt),
      individual::categorical_count_renderer_process(renderer, variables$state, categories = variables$states$get_categories())
    )

    setup_events_nimue(parameters = parameters,events = events,variables = variables,dt = dt)

    simulation_loop_nimue(
      variables = variables,
      events = events,
      processes = processes,
      timesteps = timesteps,
      progress = FALSE
    )

    df <- renderer$to_dataframe()
    df$repetition <- x
    return(df)
  })
)

# summarize safir reps
safir_reps_dt <- as.data.table(do.call(rbind,safir_reps))
safir_reps_dt[, IMild_count := IMild_count + IAsymp_count]
safir_reps_dt[, IAsymp_count := NULL]
safir_reps_dt[, IICU_count := IMVNotGetDie_count + IMVNotGetLive_count + IMVGetLive_count + IMVGetDie_count]
safir_reps_dt[, c("IMVNotGetDie_count","IMVNotGetLive_count","IMVGetLive_count","IMVGetDie_count") := NULL]
safir_reps_dt[, IHospital := IOxGetDie_count + IOxNotGetDie_count + IOxNotGetLive_count + IOxGetLive_count]
safir_reps_dt[, c("IOxGetDie_count", "IOxNotGetDie_count", "IOxNotGetLive_count", "IOxGetLive_count") := NULL]
safir_reps_dt <- melt(safir_reps_dt,id.vars = c("timestep","repetition"),variable.name = "compartment",value.name = "value")
safir_reps_dt[, compartment := gsub("(^)(\\w*)(_count)", "\\2", compartment)]
safir_reps_dt[, "t" := timestep * dt]
safir_reps_dt[, timestep := NULL]
safir_reps_dt <- safir_reps_dt[, .(ymin = quantile(value,0.025), ymax = quantile(value,0.975), y = mean(value)), by = .(t,compartment)]
safir_reps_dt[ ,"model" := "safir"]

# summarize nimue
nimue_compartments<- c("S","E","D","R","IMild","ICase","IRec","IICU","IHospital")
nimue_dt <- as.data.table(format(increasing, compartments = nimue_compartments))
nimue_dt[ , replicate := NULL]
nimue_dt <- nimue_dt[compartment %in% nimue_compartments, ]
nimue_dt[ ,"model" := "nimue"]
setnames(nimue_dt,old="value",new="y")

# combine and plot
combined_dt <- rbind(safir_reps_dt,nimue_dt,fill=TRUE)


ggplot(data = combined_dt) +
  geom_line(aes(x=t,y=y,color=compartment,linetype=model,group=model),alpha=0.85) +
  geom_ribbon(aes(x=t,ymin=ymin,ymax=ymax,fill=compartment,linetype=model),alpha=0.2) +
  facet_wrap(.~compartment,scales="free_y") +
  guides(color = FALSE,fill=FALSE)+
  theme_bw()+
  theme(strip.text.x = element_text(size=10,face = "bold"))

