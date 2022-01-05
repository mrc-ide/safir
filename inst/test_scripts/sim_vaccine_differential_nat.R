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
pop <- safir::get_population(iso3c)
pop$n <- as.integer(pop$n / 1e3)
contact_mat <- squire::get_mixing_matrix(iso3c = iso3c)

tmax <- 100
dt <- 0.5
R0 <- 4

# vaccine dosing
vaccine_doses <- 2
dose_period <- c(NaN, 28)
vaccine_set <- c(0, seq(from = 1e3, to = 1e4, length.out = tmax-1))
vaccine_set <- floor(vaccine_set)

# vaccine strategy
vaccine_coverage_mat <- nimue::strategy_matrix(strategy = "Elderly",max_coverage = 0.2)
next_dose_priority <- matrix(data = 0, nrow = vaccine_doses - 1,ncol = ncol(vaccine_coverage_mat))
next_dose_priority[1, 15:17] <- 1 # prioritize 3 oldest age groups for next dose

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
ab_parameters <- get_vaccine_ab_titre_parameters(vaccine = "Pfizer", max_dose = vaccine_doses,correlated = FALSE)

# combine parameters and verify
parameters <- make_vaccine_parameters(
  safir_parameters = parameters,
  vaccine_ab_parameters = ab_parameters,
  vaccine_set = vaccine_set,
  dose_period = dose_period,
  strategy_matrix = vaccine_coverage_mat,
  next_dose_priority_matrix = next_dose_priority
)

# ab boost for each infection
parameters$mu_ab_infection <- ab_parameters$mu_ab
parameters$mu_ab_infection <- rep(3,5)
parameters$nt_efficacy_transmission <- TRUE

parameters <- make_independent_vaccine_infection_nat_parameters(parameters = parameters, dr_vec_doses = replicate(vaccine_doses, parameters$dr_vec), dr_vec_inf = parameters$dr_vec, max_ab_inf = 8)

# create variables
timesteps <- parameters$time_period/dt
variables <- create_variables(pop = pop, parameters = parameters)
variables <- create_vaccine_variables(variables = variables,parameters = parameters)
variables <- create_natural_immunity_variables(variables = variables, parameters = parameters)
variables <- create_independent_nat_variables(variables = variables, parameters = parameters)

# make renderers
renderer <- Render$new(parameters$time_period)
dose_renderer <- Render$new(parameters$time_period)
inf_renderer <- Render$new(parameters$time_period)

nat_renderer <- Render$new(parameters$time_period)
nat_inf_renderer <- Render$new(parameters$time_period)
dose_renderer <- Render$new(parameters$time_period)

hosp_render <- create_hosp_renderers(parameters = parameters)

double_count_render_process_daily <- function(renderer, variable, dt) {
  stopifnot(inherits(variable, "DoubleVariable"))
  stopifnot(inherits(renderer, "Render"))
  function(t) {
    if ((t * dt) %% 1 == 0) {
      day <- as.integer(t * dt)
      nat <- exp(variable$get_values())
      quantiles <- quantile(x = nat, probs = c(0.025, 0.5, 0.975))
      renderer$render(name = "q025", value = quantiles[[1]], timestep = day)
      renderer$render(name = "q5", value = quantiles[[2]], timestep = day)
      renderer$render(name = "q975", value = quantiles[[3]], timestep = day)
      renderer$render(name = "mean", value = mean(x = nat), timestep = day)
    }
  }
}

# create events
events <- create_events(parameters = parameters)
events <- create_events_vaccination(events = events,parameters = parameters)
attach_event_listeners(variables = variables,events = events,parameters = parameters, dt = dt)
attach_event_listeners_vaccination(variables = variables,events = events,parameters = parameters,dt = dt)
attach_event_listeners_independent_nat(variables = variables, events = events, parameters = parameters, dt = dt)
attach_hosp_listeners(renderers = hosp_render, events = events)

# processes
processes <- list(
  independent_ab_titre_process(parameters = parameters,variables = variables,dt = dt),
  vaccination_process(parameters = parameters,variables = variables,events = events,dt = dt),
  infection_process_vaccine_cpp(parameters = parameters,variables = variables,events = events,dt = dt),
  categorical_count_renderer_process_daily(renderer = renderer,variable = variables$states,categories = variables$states$get_categories(),dt = dt),
  double_count_render_process_daily(renderer = nat_renderer, variable = variables$ab_titre,dt = dt),
  double_count_render_process_daily(renderer = nat_inf_renderer, variable = variables$ab_titre_inf,dt = dt),
  integer_count_render_process_daily(renderer = dose_renderer,variable = variables$dose_num,margin = 0:vaccine_doses,dt = dt),
  integer_count_render_process_daily(renderer = inf_renderer,variable = variables$inf_num,margin = 0:51,dt = dt)
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



# diagnostic plots ------------------------------------------------------------

ab_titre_dt <- as.data.table(nat_renderer$to_dataframe())
setnames(ab_titre_dt, "timestep", "Day")
ab_titre_dt$NAT <- "Vaccine-derived"

ab_titre_inf_dt <- as.data.table(nat_inf_renderer$to_dataframe())
setnames(ab_titre_inf_dt, "timestep", "Day")
ab_titre_inf_dt$NAT <- "Infection-derived"

ggplot(data = rbind(ab_titre_dt, ab_titre_inf_dt)) +
  geom_line(aes(x=Day,y=mean,color=NAT)) +
  geom_ribbon(aes(x=Day,ymin=q025,ymax=q975, fill=NAT),alpha=0.15) +
  theme_bw()

ggplot(data = ab_titre_dt) +
  geom_line(aes(x=Day,y=mean,color=NAT)) +
  geom_ribbon(aes(x=Day,ymin=q025,ymax=q975, fill=NAT),alpha=0.5) +
  theme_bw()

ggplot(data = ab_titre_inf_dt) +
  geom_line(aes(x=Day,y=mean,color=NAT)) +
  geom_ribbon(aes(x=Day,ymin=q025,ymax=q975, fill=NAT),alpha=0.5) +
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

# hosp/ICU
hosp_df <- process_hosp_renderers(renderers = hosp_render, parameters = parameters)
hosp_df <- melt(hosp_df, id.vars = "day")

ggplot(data = hosp_df) +
  geom_line(aes(x=day,y=value,color=variable))
