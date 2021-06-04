rm(list=ls());gc();dev.off()
library(safir)
library(individual)

N <- 1e5
dt <- 0.1

ages <- individual::IntegerVariable$new(initial_values = sample(x = 1:17,size = N,replace = T))

s0 <- rep("S",N)
s0[sample.int(n = N,size = 100,replace = TRUE)] <- sample(x = c("IMild", "IAsymp", "ICase"),size = 100,replace = T)

parameters <- safir:::get_parameters(iso3c = "ATG",time_period = 1e4)


# time step dt --------------------------------------------------------------------------------
states <- individual::CategoricalVariable$new(categories = c("S", "E", "IMild", "IAsymp", "ICase"),initial_values = s0)
times <- individual::IntegerVariable$new(initial_values = rep(0,N))

variables <- list(
  states = states,
  time = times,
  discrete_age = ages
)

events <- list(
  exposure = individual::TargetedEvent$new(N),
  mild_infection = individual::TargetedEvent$new(N),
  asymp_infection = individual::TargetedEvent$new(N),
  severe_infection = individual::TargetedEvent$new(N)
)

# Exposure events
events$exposure$add_listener(
  safir:::create_infection_update_listener(
    variables,
    "E"
  )
)

events$exposure$add_listener(
  safir:::create_exposure_scheduler_listener(
    events = events,
    variables = variables,
    parameters = parameters,
    dt = dt,
    shift = 1
  )
)

events$exposure$add_listener(
  function(timestep, target){
    times$queue_update(values = timestep,index = target)
  }
)

# IMild events
events$mild_infection$add_listener(
  safir:::create_infection_update_listener(
    variables,
    "IMild"
  )
)

events$mild_infection$add_listener(
  function(timestep, target){
    Etimes <- times$get_values(index = target)
    dwell <- timestep - Etimes
    times$queue_update(values = dwell,index = target)
  }
)


# IAsymp events
events$asymp_infection$add_listener(
  safir:::create_infection_update_listener(
    variables,
    "IAsymp"
  )
)

events$asymp_infection$add_listener(
  function(timestep, target){
    Etimes <- times$get_values(index = target)
    dwell <- timestep - Etimes
    times$queue_update(values = dwell,index = target)
  }
)

# ICase events
events$severe_infection$add_listener(
  safir:::create_infection_update_listener(
    variables,
    "ICase"
  )
)

events$severe_infection$add_listener(
  function(timestep, target){
    Etimes <- times$get_values(index = target)
    dwell <- timestep - Etimes
    times$queue_update(values = dwell,index = target)
  }
)


inf_dt <- safir:::infection_process_zzz(parameters = parameters,variables = variables,events = events,dt = dt,vaccines = NULL)

t <- 1
while (variables$states$get_size_of(c("S","E")) > 0) {
  inf_dt(timestep = t)
  for(event in events){
    event$.process()
  }
  for(var in variables){
    var$.update()
  }
  for(event in events){
    event$.tick()
  }
  t <- t + 1
  cat("t: ",t,"\n")
}

times_dt <- times$get_values()
states_dt <- sapply(X = c("IMild", "IAsymp", "ICase"),FUN = function(x){states$get_size_of(x)},USE.NAMES = TRUE)

# hist(times_dt*dt,breaks = 20)
# hist(rgamma(n = N,shape = 2,rate = 2/parameters$dur_E),breaks = 20)


# time step 1 --------------------------------------------------------------------------------
states <- individual::CategoricalVariable$new(categories = c("S", "E", "IMild", "IAsymp", "ICase"),initial_values = s0)
times <- individual::IntegerVariable$new(initial_values = rep(0,N))

variables <- list(
  states = states,
  time = times,
  discrete_age = ages
)

events <- list(
  exposure = individual::TargetedEvent$new(N),
  mild_infection = individual::TargetedEvent$new(N),
  asymp_infection = individual::TargetedEvent$new(N),
  severe_infection = individual::TargetedEvent$new(N)
)

# Exposure events
events$exposure$add_listener(
  safir:::create_infection_update_listener(
    variables,
    "E"
  )
)

events$exposure$add_listener(
  safir:::create_exposure_update_listener(
    events = events,
    variables = variables,
    parameters = parameters
  )
)

events$exposure$add_listener(
  function(timestep, target){
    times$queue_update(values = timestep,index = target)
  }
)

# IMild events
events$mild_infection$add_listener(
  safir:::create_infection_update_listener(
    variables,
    "IMild"
  )
)

events$mild_infection$add_listener(
  function(timestep, target){
    Etimes <- times$get_values(index = target)
    dwell <- timestep - Etimes
    times$queue_update(values = dwell,index = target)
  }
)


# IAsymp events
events$asymp_infection$add_listener(
  safir:::create_infection_update_listener(
    variables,
    "IAsymp"
  )
)

events$asymp_infection$add_listener(
  function(timestep, target){
    Etimes <- times$get_values(index = target)
    dwell <- timestep - Etimes
    times$queue_update(values = dwell,index = target)
  }
)

# ICase events
events$severe_infection$add_listener(
  safir:::create_infection_update_listener(
    variables,
    "ICase"
  )
)

events$severe_infection$add_listener(
  function(timestep, target){
    Etimes <- times$get_values(index = target)
    dwell <- timestep - Etimes
    times$queue_update(values = dwell,index = target)
  }
)


inf_dt <- safir:::infection_process(parameters = parameters,variables = variables,events = events)

t <- 1
while (variables$states$get_size_of(c("S","E")) > 0) {
  inf_dt(timestep = t)
  for(event in events){
    event$.process()
  }
  for(var in variables){
    var$.update()
  }
  for(event in events){
    event$.tick()
  }
  t <- t + 1
  cat("t: ",t,"\n")
}

times_1 <- times$get_values()
states_1 <- sapply(X = c("IMild", "IAsymp", "ICase"),FUN = function(x){states$get_size_of(x)},USE.NAMES = TRUE)


# compare --------------------------------------------------------------------------------

library(ggplot2)
library(data.table)


data <- data.table(
  dt = as.numeric(times_dt*dt),
  day = as.numeric(times_1),
  rgamma = rgamma(n = N,shape = 2,rate = 2/parameters$dur_E)
)
data <- melt(data,variable.name = "type")

ggplot(data = data)+
  geom_histogram(aes(x=value,fill=type,color=type),alpha=0.25,position="identity",bins=30)+
  theme_bw()

ggplot(data = data)+
  geom_histogram(aes(x=value,fill=type,color=type),position="identity",bins=30)+
  facet_wrap(.~type)+
  theme_bw()

data[type=="dt",mean(value)]
data[type=="day",mean(value)]
data[type=="rgamma",mean(value)]

data[type=="dt",var(value)]
data[type=="day",var(value)]
data[type=="rgamma",var(value)]

states_1
states_dt

