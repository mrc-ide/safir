# --------------------------------------------------------------------------------
#   test vaxx model, multi dose, no types
# --------------------------------------------------------------------------------

rm(list=ls());gc()
library(individual)
library(safir)
library(nimue)

iso3c <- "GBR"
pop <- safir:::get_population(iso3c)
pop$n <- as.integer(pop$n / 1e3)
contact_mat <- squire::get_mixing_matrix(iso3c = iso3c)

tmax <- 100
dt <- 1
R0 <- 4

# vaccine dosing
vaccine_doses <- 2
dose_period <- c(NaN, 28)
vaccine_set <- c(0, seq(from = 1e3, to = 1e4, length.out = (tmax/dt)-1))
vaccine_set <- floor(vaccine_set)
# vaccine_set <- 0*vaccine_set

# vaccine strategy
vaccine_coverage_mat <- strategy_matrix(strategy = "Elderly",max_coverage = 0.5)
next_dose_priority <- matrix(data = 0, nrow = vaccine_doses - 1,ncol = ncol(vaccine_coverage_mat))
next_dose_priority[1, 15:17] <- 1 # prioritize 3 oldest age groups for next dose
storage.mode(next_dose_priority) <- "integer"

# base parameters
parameters <- safir::get_parameters(
  population = pop$n,
  contact_matrix_set = contact_mat,
  iso3c = iso3c,
  R0 = R0,
  time_period = tmax
)
parameters$IAsymp_0[1:17] <- 3 # some additional infectives at t=0
parameters$population <- parameters$population + 3
pop$n <- parameters$population # fix this. we don't want to use this twice.

# attach vaccine parameters (roll into function eventually)
parameters$N_prioritisation_steps <- nrow(vaccine_coverage_mat)
parameters$vaccine_coverage_mat <- vaccine_coverage_mat
parameters$next_dose_priority <- next_dose_priority

parameters$vaccine_set <- vaccine_set
parameters$dose_period <- dose_period
parameters$N_phase <- vaccine_doses

# attach Ab dynamics parameters
ab_parameters <- get_vaccine_ab_titre_parameters(vaccine = "Pfizer", max_dose = vaccine_doses)
parameters <- c(parameters, ab_parameters)

# create variables
timesteps <- parameters$time_period/dt
variables <- create_variables(pop = pop, parameters = parameters)
variables <- create_vaccine_variables(variables = variables,pop = parameters$population,max_dose = vaccine_doses)

# create events
events <- create_events(parameters = parameters)
events <- create_events_vaccination(events = events,parameters = parameters)
attach_event_listeners(variables = variables,events = events,parameters = parameters, dt = dt)
attach_event_listeners_vaccination(variables = variables,events = events,parameters = parameters,dt = dt)

# make renderers
renderer <- Render$new(parameters$time_period)
ef_inf_renderer <- matrix(data = NaN,nrow = parameters$time_period,ncol = sum(parameters$population))
dose_renderer <- Render$new(parameters$time_period)

double_count_render_process_daily <- function(variable, dt) {
  stopifnot(inherits(variable, "DoubleVariable"))
  function(t) {
    if ((t * dt) %% 1 == 0) {
      ef_inf_renderer[as.integer(t * dt), ] <<- variable$get_values()
    }
  }
}

# processes
processes <- list(
  vaccine_ab_titre_process(parameters = parameters,variables = variables,events = events,dt = dt),
  vaccination_process(parameters = parameters,variables = variables,events = events,dt = dt),
  infection_process_vaccine(parameters = parameters,variables = variables,events = events,dt = dt),
  categorical_count_renderer_process_daily(renderer = renderer,variable = variables$states,categories = variables$states$get_categories(),dt = dt),
  double_count_render_process_daily(variable = variables$ef_infection,dt = dt),
  integer_count_render_process_daily(renderer = dose_renderer,variable = variables$dose_num,margin = 0:2,dt = dt)
)

setup_events_vaccine(parameters = parameters,events = events,variables = variables,dt = dt)

system.time(simulation_loop_vaccine(
  variables = variables,
  events = events,
  processes = processes,
  timesteps = timesteps,
  TRUE
))




# plot
vaccinated <- variables$dose_num$get_index_of(set = 0)
vaccinated <- vaccinated$not()

ef_infection <- ef_inf_renderer[, vaccinated$to_vector()]
start <- apply(ef_infection, 2, function(x){ which(x < 1)[1] })

ef_infection <- lapply(X = 1:ncol(ef_infection),FUN = function(x){
  ef_infection[start[x]:nrow(ef_infection) , x]
})

ef_infection_dt <- data.table(
  id = rep(1:length(ef_infection), times = sapply(ef_infection,length)),
  t = unlist(lapply(ef_infection,function(x){1:length(x)})),
  ef_inf = unlist(ef_infection)
)
ef_infection_dt[, ef_inf := 1 - ef_inf]

eff_sum_dt <- ef_infection_dt[, .(mean = mean(ef_inf), lo = quantile(ef_inf,probs =0.025) , hi = quantile(ef_inf,probs = 0.975)) , by = .(t)]

ggplot(data = eff_sum_dt) +
  geom_line(aes(x=t,y=mean)) +
  geom_ribbon(aes(x=t,ymin=lo,ymax=hi),alpha=0.5) +
  theme_bw()



dose_out <- dose_renderer$to_dataframe()
matplot(dose_out[,-1],type="l",lty=1)


saf_dt <- as.data.table(renderer$to_dataframe())
saf_dt[, IMild_count := IMild_count + IAsymp_count]
saf_dt[, IAsymp_count := NULL]
saf_dt <- melt(saf_dt,id.vars = c("timestep"),variable.name = "name")
saf_dt[, name := gsub("(^)(\\w*)(_count)", "\\2", name)]
setnames(x = saf_dt,old = c("timestep","name","value"),new = c("t","compartment","y"))

ggplot(data = saf_dt, aes(t,y,color = compartment)) +
  geom_line() +
  geom_line() +
  facet_wrap(~compartment, scales = "free")

