rm(list=ls());gc()

# 1. give out n+1 prioritized doses before dose n
# 2. make stepping to phase n+1 depend on prioritized coverage too

vaccine_doses <- 4

max_coverage_d2 = 0.9
max_coverage_d3 = 0.8
max_coverage_d4 = 0.5

age_groups_covered_d2 = 15
age_groups_covered_d3 = 14
age_groups_covered_d4 = 5

vaccine_coverage_strategy <- list()
vaccine_coverage_strategy[[1]] <-
  nimue::strategy_matrix(strategy = "Elderly", max_coverage = max_coverage_d2)[1:age_groups_covered_d2, ]

vaccine_coverage_strategy[[2]] <-
  nimue::strategy_matrix(strategy = "Elderly", max_coverage = max_coverage_d2)[1:age_groups_covered_d2, ]

vaccine_coverage_strategy[[3]] <-
  nimue::strategy_matrix(strategy = "Elderly", max_coverage = max_coverage_d3)[1:age_groups_covered_d3, ]

vaccine_coverage_strategy[[4]] <-
  nimue::strategy_matrix(strategy = "Elderly", max_coverage = max_coverage_d4)[1:age_groups_covered_d4, ]


# reduce coverage in some younger age groups
vaccine_coverage_strategy[[1]][,c((17 - age_groups_covered_d2 + 1):5)] <- max_coverage_d2 * .83 * vaccine_coverage_strategy[[1]][,c((17 - age_groups_covered_d2 + 1):5)]
vaccine_coverage_strategy[[2]][,c((17 - age_groups_covered_d2 + 1):5)] <- max_coverage_d2 * .83 * vaccine_coverage_strategy[[2]][,c((17 - age_groups_covered_d2 + 1):5)]
vaccine_coverage_strategy[[3]][,c((17 - age_groups_covered_d3 + 1):10)] <- max_coverage_d3 * .625 * vaccine_coverage_strategy[[3]][,c((17 - age_groups_covered_d3 + 1):10)]

next_dose_priority <- matrix(data = 0, nrow = vaccine_doses - 1, ncol = ncol(vaccine_coverage_strategy[[1]]))

# next_dose_priority[1,(17 - age_groups_covered_d2 + 1):17] <- 1
# next_dose_priority[2,(17 - age_groups_covered_d3 + 1):17] <- 1
# next_dose_priority[3,(17 - age_groups_covered_d4 + 1):17] <- 1

library(individual)
library(data.table)
library(ggplot2)
library(safir)
library(nimue)

iso3c <- "GBR"
pop <- safir:::get_population(iso3c)
pop$n <- as.integer(pop$n / 5e2)
contact_mat <- squire::get_mixing_matrix(iso3c = iso3c)

tmax <- 400
dt <- 0.5
R0 <- 4

# vaccine dosing
dose_period <- c(NaN, 10, 10, 10)
vaccine_set <- c(0, rep(1e3, tmax - 1))

parameters <- safir::get_parameters(
  population = pop$n,
  contact_matrix_set = contact_mat,
  iso3c = iso3c,
  R0 = R0,
  time_period = tmax,
  dt = dt
)

mu_ab_list_4dose <- data.frame(
  name = c("Pfizer", "AstraZeneca", "Sinovac", "Moderna"),
  d1 = c(13/94, 1/59, 28/164, ((185 + 273)/2)/321),
  d2 = c(223/94, 32/59,28/164, 654/158)
)
mu_ab_list_4dose <- cbind(mu_ab_list_4dose, d3 = mu_ab_list_4dose[, 3], d4 = mu_ab_list_4dose[, 3])

ab_parameters <- get_vaccine_ab_titre_parameters(vaccine = "Pfizer", max_dose = vaccine_doses,correlated = FALSE, mu_ab_list = mu_ab_list_4dose)

# combine parameters and verify
parameters <- make_vaccine_parameters(
  safir_parameters = parameters,
  vaccine_ab_parameters = ab_parameters,
  vaccine_set = vaccine_set,
  dose_period = dose_period,
  strategy_matrix = vaccine_coverage_strategy,
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
  integer_count_render_process_daily(renderer = dose_renderer,variable = variables$dose_num,margin = 0:vaccine_doses,dt = dt)
)

setup_events(parameters = parameters,events = events,variables = variables,dt = dt)

doses_tracker <- cbind("dose1"=rep(NaN,timesteps),"dose2"=rep(NaN,timesteps),"dose3"=rep(NaN,timesteps),"dose4"=rep(NaN,timesteps))
doses_tracker <- as.data.frame(doses_tracker)

# for that to work, surround any call to `assign_doses` with:

# doses_pre <- doses_left # EXTRA
#
# assign_doses
#
# doses_post <- doses_left # EXTRA
# doses_given <- doses_pre - doses_post # EXTRA
# .GlobalEnv$doses_tracker[day, phase+1] <- doses_given # EXTRA

system.time(simulation_loop_safir(
  variables = variables,
  events = events,
  processes = processes,
  timesteps = timesteps,
  variables_dont_update = c("discrete_age", "phase"),
  progress = TRUE
))


# plot: vaccinations
dose_out <- dose_renderer$to_dataframe()
colnames(dose_out)[2:(vaccine_doses+2)] <- as.character(0:vaccine_doses)
dose_out <- melt(as.data.table(dose_out),id.vars="timestep")
setnames(dose_out, "variable", "dose")

ggplot(data = dose_out) +
  geom_line(aes(x=timestep,y=value,color=dose)) +
  xlab("day") +
  theme_bw()


# plot cumulative vaxx


doses_cum <- as.data.table(doses_tracker)
for (i in names(doses_cum)) {
  doses_cum[is.nan(get(i)), (i):=0]
}
doses_cum[, "dose1" := cumsum(dose1)]
doses_cum[, "dose2" := cumsum(dose2)]
doses_cum[, "dose3" := cumsum(dose3)]
doses_cum[, "dose4" := cumsum(dose4)]
doses_cum[, "day" := 1:nrow(doses_cum)]

doses_cum <- melt(doses_cum, id.vars="day")

ggplot(data = doses_cum, aes(x = day, y = value, col = variable)) +
  geom_line(size = 0.8) +
  #geom_bar(stat = "identity") +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.text.align = 0) +
  scale_color_viridis_d(option = "C", end = 0.9) +
  labs(x = "timestep", y = "cumulative doses", col = "dose number")
