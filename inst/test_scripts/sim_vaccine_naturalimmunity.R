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
vaccine_set <- vaccine_set*0

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
inf_renderer <- Render$new(parameters$time_period)

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
  integer_count_render_process_daily(renderer = dose_renderer,variable = variables$dose_num,margin = 0:vaccine_doses,dt = dt),
  integer_count_render_process_daily(renderer = inf_renderer,variable = variables$inf_num,margin = 0:51,dt = dt)
)

setup_events_vaccine(parameters = parameters,events = events,variables = variables,dt = dt)

# # stuff we'd like to check
# debug(processes[[1]])
# debugonce(events$recovery$.listeners[[2]])
# debug(events$recovery$.listeners[[3]])

system.time(simulation_loop_safir(
  variables = variables,
  events = events,
  processes = processes,
  timesteps = timesteps,
  variables_dont_update = c("discrete_age", "phase"),
  progress = TRUE
))



# diagnostic plots ------------------------------------------------------------

# plot: ab titre
vaccinated_or_infected <- variables$dose_num$get_index_of(set = 0)
vaccinated_or_infected$not(inplace = TRUE)

vaccinated_or_infected$or(variables$inf_num$get_index_of(set = 0)$not(inplace = TRUE))

ab_titre <- ab_renderer[, vaccinated_or_infected$to_vector()]
ab_titre[which(!is.finite(ab_titre))] <- NaN
start <- apply(ab_titre, 2, function(x){ which(abs(x - 0) > 2e-7)[1] })

ab_titre <- lapply(X = 1:ncol(ab_titre),FUN = function(x){
  ab_titre[start[x]:nrow(ab_titre) , x]
})

ab_titre_dt <- data.table(
  ab = unlist(ab_titre),
  t = unlist(lapply(ab_titre, function(ab_titre){1:length(ab_titre)})),
  id = rep(1:length(ab_titre), times = vapply(ab_titre,length,integer(1)))
)

ab_titre_quant_dt <- ab_titre_dt[, .(mean = mean(ab), lo = quantile(ab,probs =0.025) , hi = quantile(ab,probs = 0.975)) , by = .(t)]

ggplot(data = ab_titre_quant_dt) +
  geom_line(aes(x=t,y=mean)) +
  geom_ribbon(aes(x=t,ymin=lo,ymax=hi),alpha=0.5) +
  theme_bw()

# plot: vaccinations
dose_out <- dose_renderer$to_dataframe()
colnames(dose_out)[2:(vaccine_doses+2)] <- as.character(0:vaccine_doses)
dose_out <- as.data.table(dose_out)
dose_out <- melt(dose_out, id.vars="timestep")
setnames(dose_out, "variable", "dose")

ggplot(data = dose_out) +
  geom_line(aes(x=timestep,y=value,color=dose)) +
  theme_bw()

# states
saf_dt <- as.data.table(renderer$to_dataframe())
saf_dt[, IMild_count := IMild_count + IAsymp_count]
saf_dt[, IAsymp_count := NULL]
saf_dt <- melt(saf_dt,id.vars = c("timestep"),variable.name = "name")
saf_dt[, name := gsub("(^)(\\w*)(_count)", "\\2", name)]
setnames(x = saf_dt,old = c("timestep","name","value"),new = c("t","compartment","y"))

ggplot(data = saf_dt, aes(t,y,color = compartment)) +
  geom_line() +
  facet_wrap(~compartment, scales = "free")

# infections
inf_out <- inf_renderer$to_dataframe()
inf_out <- inf_out[, -which(sapply(inf_out, function(x){all(x == 0)}, USE.NAMES = FALSE))]
inf_out <- as.data.table(inf_out)
inf_out <- melt(inf_out, id.vars = c("timestep"),variable.name = "infections")
inf_out[, infections := gsub("(^X)(\\w*)(_count)", "\\2", infections)]
inf_out[, t := timestep * dt]
inf_out[, timestep := NULL]

ggplot(data = inf_out, aes(t,value,color = infections)) +
  geom_line(size = 1.25, alpha = 0.75)



# gut check

not_inf_or_vaccinated <- variables$inf_num$get_index_of(set = 0)$and(variables$dose_num$get_index_of(set = 0))
all(!is.finite(variables$ab_titre$get_values(not_inf_or_vaccinated)))
