# --------------------------------------------------------------------------------
#   test vaxx model, multi dose, no types
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
vaccine_doses <- 3
dose_period <- c(NaN, 28, 56)
# vaccine_set <- vaccine_set <- c(0, seq(from = 1e3, to = 1e4, length.out = tmax-1))
# vaccine_set <- floor(vaccine_set)
vaccine_set <- c(0, rep(1e3, tmax - 1))

# vaccine strategy
vaccine_coverage_mat <- strategy_matrix(strategy = "All",max_coverage = 0.5)
next_dose_priority <- matrix(data = 0, nrow = vaccine_doses - 1,ncol = ncol(vaccine_coverage_mat))
next_dose_priority[1, 15:17] <- 1 # prioritize 3 oldest age groups for next dose
next_dose_priority[2, 15:17] <- 1 # prioritize 3 oldest age groups for next dose

# base parameters
parameters <- safir::get_parameters(
  population = pop$n,
  contact_matrix_set = contact_mat,
  iso3c = iso3c,
  R0 = R0,
  time_period = tmax,
  dt = dt
)

# vaccine parameters
mu_ab_list_3dose <- data.frame(
  name = c("Pfizer", "AstraZeneca", "Sinovac", "Moderna"),
  d1 = c(13/94, 1/59, 28/164, ((185 + 273)/2)/321),
  d2 = c(223/94, 32/59,28/164, 654/158)
)
mu_ab_list_3dose <- cbind(mu_ab_list_3dose, d3 = mu_ab_list_3dose[, 3])

ab_parameters <- get_vaccine_ab_titre_parameters(vaccine = "Pfizer", max_dose = vaccine_doses,correlated = TRUE,mu_ab_list = mu_ab_list_3dose)

# combine parameters and verify
parameters <- make_vaccine_parameters(
  safir_parameters = parameters,
  vaccine_ab_parameters = ab_parameters,
  vaccine_set = vaccine_set,
  dose_period = dose_period,
  strategy_matrix = vaccine_coverage_mat,
  next_dose_priority_matrix = next_dose_priority
)

# create variables
timesteps <- parameters$time_period/dt
variables <- create_variables(pop = pop, parameters = parameters)
variables <- create_vaccine_variables(variables = variables,parameters = parameters)

# create events
events <- create_events(parameters = parameters)
events <- create_events_vaccination(events = events,parameters = parameters)
attach_event_listeners(variables = variables,events = events,parameters = parameters, dt = dt)
attach_event_listeners_vaccination(variables = variables,events = events,parameters = parameters,dt = dt)

# make renderers
renderer <- Render$new(parameters$time_period)
nat_renderer <- Render$new(parameters$time_period)
dose_renderer <- Render$new(parameters$time_period)

double_count_render_process_daily <- function(renderer, variable, dt) {
  stopifnot(inherits(variable, "DoubleVariable"))
  stopifnot(inherits(renderer, "Render"))
  function(t) {
    if ((t * dt) %% 1 == 0) {
      day <- as.integer(t * dt)
      nat <- exp(variable$get_values())
      renderer$render(name = "q025", value = quantile(x = nat, probs = 0.025), timestep = day)
      renderer$render(name = "q975", value = quantile(x = nat, probs = 0.975), timestep = day)
      renderer$render(name = "q5", value = quantile(x = nat, probs = 0.5), timestep = day)
      renderer$render(name = "mean", value = mean(x = nat), timestep = day)
    }
  }
}

# processes
processes <- list(
  vaccine_ab_titre_process(parameters = parameters,variables = variables,dt = dt),
  vaccination_process(parameters = parameters,variables = variables,events = events,dt = dt),
  infection_process_vaccine_cpp(parameters = parameters,variables = variables,events = events,dt = dt),
  categorical_count_renderer_process_daily(renderer = renderer,variable = variables$states,categories = variables$states$get_categories(),dt = dt),
  double_count_render_process_daily(renderer = nat_renderer, variable = variables$ab_titre, dt = dt),
  integer_count_render_process_daily(renderer = dose_renderer,variable = variables$dose_num,margin = 0:vaccine_doses,dt = dt)
)

setup_events(parameters = parameters,events = events,variables = variables,dt = dt)

system.time(simulation_loop_safir(
  variables = variables,
  events = events,
  processes = processes,
  timesteps = timesteps,
  variables_dont_update = c("discrete_age", "phase"),
  progress = TRUE
))


# # plot: ab titre
# ab_titre_dt <- as.data.table(nat_renderer$to_dataframe())
# setnames(ab_titre_dt, "timestep", "Day")
#
# ggplot(data = ab_titre_dt) +
#   geom_line(aes(x=Day,y=mean)) +
#   geom_ribbon(aes(x=Day,ymin=q025,ymax=q975),alpha=0.5) +
#   theme_bw()

# plot: vaccinations
dose_out <- dose_renderer$to_dataframe()
colnames(dose_out)[2:(vaccine_doses+2)] <- as.character(0:vaccine_doses)
dose_out <- melt(as.data.table(dose_out),id.vars="timestep")
setnames(dose_out, "variable", "dose")

ggplot(data = dose_out) +
  geom_line(aes(x=timestep,y=value,color=dose)) +
  theme_bw()

# # cumulative vaccinations
# dcast(dose_out, timestep ~ dose, value.var="value", fill=0L)
#
# dose_cum <- dose_out[, "cumulative" := cumsum(value), by = "dose"]
#
#
# ggplot(data = dose_cum[dose!=0,]) +
#   geom_line(aes(x=timestep,y=cumulative,color=dose)) +
#   theme_bw()


# # states
# saf_dt <- as.data.table(renderer$to_dataframe())
# saf_dt[, IMild_count := IMild_count + IAsymp_count]
# saf_dt[, IAsymp_count := NULL]
# saf_dt <- melt(saf_dt,id.vars = c("timestep"),variable.name = "name")
# saf_dt[, name := gsub("(^)(\\w*)(_count)", "\\2", name)]
# setnames(x = saf_dt,old = c("timestep","name","value"),new = c("t","compartment","y"))
#
# ggplot(data = saf_dt, aes(t,y,color = compartment)) +
#   geom_line() +
#   geom_line() +
#   facet_wrap(~compartment, scales = "free")

