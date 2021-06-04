rm(list=ls());gc();dev.off()
library(safir)
library(individual)

N <- 1e4

ages_1 <- individual::IntegerVariable$new(initial_values = sample(x = 1:17,size = N,replace = T))
ages_dt <- individual::IntegerVariable$new(initial_values = ages_1$get_values())

s0 <- rep("S",N)
s0[sample.int(n = N,size = 100,replace = TRUE)] <- sample(x = c("IMild", "IAsymp", "ICase"),size = 100,replace = T)

states_1 <- individual::CategoricalVariable$new(categories = c("S", "E", "IMild", "IAsymp", "ICase"),initial_values = s0)
states_dt <- individual::CategoricalVariable$new(categories = c("S", "E", "IMild", "IAsymp", "ICase"),initial_values = s0)

variables_1 <- list(states=states_1, discrete_age = ages_1)
variables_dt <- list(states=states_dt, discrete_age = ages_dt)

exposure_1 = individual::TargetedEvent$new(N)
exposure_dt = individual::TargetedEvent$new(N)

events_1 <- list(exposure = exposure_1)
events_dt <- list(exposure = exposure_dt)

steps_1 <- individual::IntegerVariable$new(initial_values = rep(0,N))
steps_dt <- individual::IntegerVariable$new(initial_values = rep(0,N))

parameters <- safir:::get_parameters(iso3c = "ATG",time_period = 1e4)
parameters$beta_set <- parameters$beta_set * 2

exposure_1$add_listener(
  safir:::create_infection_update_listener(
    variables_1,
    "E"
  )
)

exposure_1$add_listener(
  function(t, target){
    steps_1$queue_update(values = t,index = target)
  }
)

exposure_dt$add_listener(
  safir:::create_infection_update_listener(
    variables_dt,
    "E"
  )
)

exposure_dt$add_listener(
  function(t, target){
    steps_dt$queue_update(values = t,index = target)
  }
)


# daily time step
inf_1 <- safir:::infection_process(parameters = parameters,variables = variables_1,events = events_1)

t <- 1
while (variables_1$states$get_size_of('S') > 0) {
  if(t==10002){
    browser()
  }
  inf_1(timestep = t)
  exposure_1$.process()
  steps_1$.update()
  variables_1$states$.update()
  exposure_1$.tick()
  t <- t + 1
  cat("t: ",t,"\n")
}

hist(steps_1$get_values())


# user select time step
dt <- 0.1

inf_dt <- safir:::infection_process_zzz(parameters = parameters,variables = variables_dt,events = events_dt,dt = dt,vaccines = NULL)

t <- 1
while (variables_dt$states$get_size_of('S') > 0) {
  # if(t==10002){
  #   browser()
  # }
  inf_dt(timestep = t)
  exposure_dt$.process()
  steps_dt$.update()
  variables_dt$states$.update()
  exposure_dt$.tick()
  t <- t + 1
  cat("t: ",t,"\n")
}

hist(steps_dt$get_values() * dt)
