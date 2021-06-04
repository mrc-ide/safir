rm(list=ls());gc()
library(safir)
library(individual)

N <- 5e5
dt <- 0.1

ages <- individual::IntegerVariable$new(initial_values = sample(x = 1:17,size = N,replace = T))

s0 <- rep("S",N)
s0[sample.int(n = N,size = 100,replace = TRUE)] <- sample(x = c("IMild", "IAsymp", "ICase"),size = 100,replace = T)

parameters <- safir:::get_parameters(iso3c = "ATG",time_period = 1e4)
parameters$beta_set <- parameters$beta_set*1e2

exp_shift<-1
other_shift<-1

exp_func <- safir:::make_rexp_simple
# exp_func <- safir:::make_rexp

# ----------------------------------------------------------------------------------------------------
# time step dt --------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------

states <- individual::CategoricalVariable$new(categories = c("S", "E", "IMild", "IAsymp", "ICase","R"),initial_values = s0)
E_times <- individual::IntegerVariable$new(initial_values = rep(-1,N))
I_MildAsymp_times <- individual::IntegerVariable$new(initial_values = rep(-1,N))
I_severe_times <- individual::IntegerVariable$new(initial_values = rep(-1,N))

variables <- list(
  states = states,
  E_time = E_times,
  I_MildAsymp = I_MildAsymp_times,
  I_severe = I_severe_times,
  discrete_age = ages
)

events <- list(
  exposure = individual::TargetedEvent$new(N),
  mild_infection = individual::TargetedEvent$new(N),
  asymp_infection = individual::TargetedEvent$new(N),
  severe_infection = individual::TargetedEvent$new(N),
  recovery_from_mild_asymp = individual::TargetedEvent$new(N),
  recovery_from_severe = individual::TargetedEvent$new(N)
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
    shift = exp_shift
  )
)

events$exposure$add_listener(
  function(timestep, target){
    E_times$queue_update(values = timestep,index = target)
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
    Etimes <- E_times$get_values(index = target)
    dwell <- timestep - Etimes
    E_times$queue_update(values = dwell,index = target)

    I_MildAsymp_times$queue_update(values = timestep,index = target)
  }
)

events$mild_infection$add_listener(
  safir:::create_event_scheduler_listener(
    event = events$recovery_from_mild_asymp,
    duration = parameters$dur_IMild,
    func = exp_func,
    shift = other_shift,
    dt = dt
  )
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
    Etimes <- E_times$get_values(index = target)
    dwell <- timestep - Etimes
    E_times$queue_update(values = dwell,index = target)

    I_MildAsymp_times$queue_update(values = timestep,index = target)
  }
)

events$asymp_infection$add_listener(
  safir:::create_event_scheduler_listener(
    event = events$recovery_from_mild_asymp,
    duration = parameters$dur_IAsymp,
    func = exp_func,
    shift = other_shift,
    dt = dt
  )
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
    Etimes <- E_times$get_values(index = target)
    dwell <- timestep - Etimes
    E_times$queue_update(values = dwell,index = target)

    I_severe_times$queue_update(values = timestep,index = target)
  }
)

events$severe_infection$add_listener(
  safir:::create_event_scheduler_listener(
    event = events$recovery_from_severe,
    duration = parameters$dur_ICase,
    func = safir:::make_rerlang,
    shift = other_shift,
    dt = dt
  )
)


# recovery from mild/asymp
events$recovery_from_mild_asymp$add_listener(
  safir:::create_infection_update_listener(
    variables,
    "R"
  )
)

events$recovery_from_mild_asymp$add_listener(
  function(timestep, target){
    t0 <- I_MildAsymp_times$get_values(index = target)
    dwell <- timestep - t0
    I_MildAsymp_times$queue_update(values = dwell,index = target)
  }
)

# recovery from severe
events$recovery_from_severe$add_listener(
  safir:::create_infection_update_listener(
    variables,
    "R"
  )
)

events$recovery_from_severe$add_listener(
  function(timestep, target){
    t0 <- I_severe_times$get_values(index = target)
    dwell <- timestep - t0
    I_severe_times$queue_update(values = dwell,index = target)
  }
)

#  --------------------------------------------------------------------------------
# initialize disease progression

# 1
target_mild <- states$get_index_of("IMild")

init_mild <- safir:::create_event_scheduler_listener(
  event = events$recovery_from_mild_asymp,
  duration = parameters$dur_IMild,
  func = exp_func,
  shift = other_shift,
  dt = dt
)

I_MildAsymp_times$queue_update(values = 1,index = target_mild)

init_mild(timestep = 1,target = target_mild)

# 2
target_asymp <- states$get_index_of("IAsymp")

init_asymp <- safir:::create_event_scheduler_listener(
  event = events$recovery_from_mild_asymp,
  duration = parameters$dur_IAsymp,
  func = exp_func,
  shift = other_shift,
  dt = dt
)

I_MildAsymp_times$queue_update(values = 1,index = target_asymp)

init_asymp(timestep = 1,target = target_asymp)

# 3
target_case <- states$get_index_of("ICase")

init_case <- safir:::create_event_scheduler_listener(
  event = events$recovery_from_severe,
  duration = parameters$dur_ICase,
  func = exp_func,
  shift = other_shift,
  dt = dt
)

I_severe_times$queue_update(values = 1,index = target_case)

init_case(timestep = 1,target = target_case)

# sim ----------------------------------------------------------------------------------------

inf_dt <- safir:::infection_process_zzz(parameters = parameters,variables = variables,events = events,dt = dt,vaccines = NULL)

# lapply(events$mild_infection$.listeners,debugonce)

t <- 1
# while (variables$states$get_size_of(c("S","E")) > 0) {
while (variables$states$get_size_of(c("S","E", "ICase", "IMild", "IAsymp")) > 0) {
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
  # cat("t: ",t," size of states: ",variables$states$get_size_of(c("S","E","IMild", "IAsymp", "ICase")) ,"\n")
  cat("t: ",t," size of states: ",variables$states$get_size_of(c("S","E")) ,"\n")
}

states_dt <- sapply(X = states$get_categories(),FUN = function(x){states$get_size_of(x)},USE.NAMES = TRUE)

E_times_dt <- E_times$get_values()
E_times_dt <- E_times_dt[E_times_dt >= 0]

I_severe_times_dt <- I_severe_times$get_values()
I_severe_times_dt <- I_severe_times_dt[I_severe_times_dt >=0]

I_mild_times_dt <- I_MildAsymp_times$get_values()
I_mild_times_dt <- I_mild_times_dt[I_mild_times_dt>=0]

# par(mfrow=c(3,1))
# hist(E_times_dt*dt,breaks = 30,main=paste0("dt: ",dt," historgram of E times, mean: ",mean(E_times_dt*dt)," expected: ",parameters$dur_E))
# hist(I_severe_times_dt*dt,breaks=30,main=paste0("dt: ",dt," historgram of I severe times, mean: ",mean(I_severe_times_dt*dt)," expected: ",parameters$dur_ICase))
# hist(I_mild_times_dt*dt,breaks=30,main=paste0("dt: ",dt," historgram of I mild/asymp times, mean: ",mean(I_mild_times_dt*dt)," expected: ",parameters$dur_IMild))
# par(mfrow=c(1,1))




# ----------------------------------------------------------------------------------------------------
# time step 1 --------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------

states <- individual::CategoricalVariable$new(categories = c("S", "E", "IMild", "IAsymp", "ICase","R"),initial_values = s0)
E_times <- individual::IntegerVariable$new(initial_values = rep(-1,N))
I_MildAsymp_times <- individual::IntegerVariable$new(initial_values = rep(-1,N))
I_severe_times <- individual::IntegerVariable$new(initial_values = rep(-1,N))

variables <- list(
  states = states,
  E_time = E_times,
  I_MildAsymp = I_MildAsymp_times,
  I_severe = I_severe_times,
  discrete_age = ages
)

events <- list(
  exposure = individual::TargetedEvent$new(N),
  mild_infection = individual::TargetedEvent$new(N),
  asymp_infection = individual::TargetedEvent$new(N),
  severe_infection = individual::TargetedEvent$new(N),
  recovery_from_mild_asymp = individual::TargetedEvent$new(N),
  recovery_from_severe = individual::TargetedEvent$new(N)
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
    E_times$queue_update(values = timestep,index = target)
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
    Etimes <- E_times$get_values(index = target)
    dwell <- timestep - Etimes
    E_times$queue_update(values = dwell,index = target)

    I_MildAsymp_times$queue_update(values = timestep,index = target)
  }
)

events$mild_infection$add_listener(
  safir:::create_progression_listener(
    event = events$recovery_from_mild_asymp,
    duration = parameters$dur_IMild,
    func = safir:::r_exp
  )
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
    Etimes <- E_times$get_values(index = target)
    dwell <- timestep - Etimes
    E_times$queue_update(values = dwell,index = target)

    I_MildAsymp_times$queue_update(values = timestep,index = target)
  }
)

events$asymp_infection$add_listener(
  safir:::create_progression_listener(
    event = events$recovery_from_mild_asymp,
    duration = parameters$dur_IAsymp,
    func = safir:::r_exp
  )
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
    Etimes <- E_times$get_values(index = target)
    dwell <- timestep - Etimes
    E_times$queue_update(values = dwell,index = target)

    I_severe_times$queue_update(values = timestep,index = target)
  }
)

events$severe_infection$add_listener(
  safir:::create_progression_listener(
    event = events$recovery_from_severe,
    duration = parameters$dur_ICase
  )
)

# recovery from mild/asymp
events$recovery_from_mild_asymp$add_listener(
  safir:::create_infection_update_listener(
    variables,
    "R"
  )
)

events$recovery_from_mild_asymp$add_listener(
  function(timestep, target){
    t0 <- I_MildAsymp_times$get_values(index = target)
    dwell <- timestep - t0
    I_MildAsymp_times$queue_update(values = dwell,index = target)
  }
)

# recovery from severe
events$recovery_from_severe$add_listener(
  safir:::create_infection_update_listener(
    variables,
    "R"
  )
)

events$recovery_from_severe$add_listener(
  function(timestep, target){
    t0 <- I_severe_times$get_values(index = target)
    dwell <- timestep - t0
    I_severe_times$queue_update(values = dwell,index = target)
  }
)

#  --------------------------------------------------------------------------------
# initialize disease progression

# 1
target_mild <- states$get_index_of("IMild")

init_mild <- safir:::create_progression_listener(
  event = events$recovery_from_mild_asymp,
  duration = parameters$dur_IMild,
  func = safir:::r_exp
)

I_MildAsymp_times$queue_update(values = 1,index = target_mild)

init_mild(timestep = 1,target = target_mild)

# 2
target_asymp <- states$get_index_of("IAsymp")

init_asymp <- safir:::create_progression_listener(
  event = events$recovery_from_mild_asymp,
  duration = parameters$dur_IAsymp,
  func = safir:::r_exp
)

I_MildAsymp_times$queue_update(values = 1,index = target_asymp)

init_asymp(timestep = 1,target = target_asymp)

# 3
target_case <- states$get_index_of("ICase")

init_case <- safir:::create_progression_listener(
  event = events$recovery_from_severe,
  duration = parameters$dur_ICase
)

I_severe_times$queue_update(values = 1,index = target_case)

init_case(timestep = 1,target = target_case)

#  --------------------------------------------------------------------------------
# sim

inf_1 <- safir:::infection_process(parameters = parameters,variables = variables,events = events)

t <- 1
# while (variables$states$get_size_of(c("S","E")) > 0) {
while (variables$states$get_size_of(c("S","E", "ICase", "IMild", "IAsymp")) > 0) {
  inf_1(timestep = t)
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
  cat("t: ",t," size of states: ",variables$states$get_size_of(c("S","E")) ,"\n")
}

states_1 <- sapply(X = states$get_categories(),FUN = function(x){states$get_size_of(x)},USE.NAMES = TRUE)

E_times_1 <- E_times$get_values()
E_times_1 <- E_times_1[E_times_1 >= 0]

I_severe_times_1 <- I_severe_times$get_values()
I_severe_times_1 <- I_severe_times_1[I_severe_times_1 >=0]

I_mild_times_1 <- I_MildAsymp_times$get_values()
I_mild_times_1 <- I_mild_times_1[I_mild_times_1>=0]


# compare --------------------------------------------------------------------------------

library(ggplot2)
library(data.table)

data_E <- data.table(
  dt = as.numeric(E_times_dt*dt),
  day = as.numeric(E_times_1)
)
data_E <- melt(data_E,variable.name = "timestep")
data_E[, "state" := "E"]

data_severe <- data.table(
  timestep = c(rep("dt",length(I_severe_times_dt)),rep("day",length(I_severe_times_1))),
  value = c(as.numeric(I_severe_times_dt*dt),as.numeric(I_severe_times_1))
)
data_severe[, "state" := "Severe"]

data_mild <- data.table(
  timestep = c(rep("dt",length(I_mild_times_dt)),rep("day",length(I_mild_times_1))),
  value = c(as.numeric(I_mild_times_dt*dt),as.numeric(I_mild_times_1))
)
data_mild[, "state" := "Mild/Asymp"]

data <- rbind(data_E,data_severe,data_mild)

ggplot(data = data) +
  geom_histogram(aes(x=value,fill=timestep,color=timestep,y=..density..),alpha=0.25,position="identity",bins=30)+
  facet_wrap(.~state,scales="free") +
  theme_bw()

ggplot(data = data) +
  geom_histogram(aes(x=value,fill=timestep,color=timestep,y=..density..),alpha=0.25,position="identity",bins=35)+
  facet_grid(timestep~state,scales="free") +
  theme_bw()

