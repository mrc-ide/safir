rm(list=ls());gc()
library(safir)
library(individual)
library(data.table)
library(parallel)
library(ggplot2)

iso3c <- "ATG"
pop <- safir:::get_population(iso3c)

ncores <- detectCores()-6

nrep <- 25
# Scale it for speed
# pop$n <- as.integer(pop$n / 10)

# Create our simulation parameters
R0 <- 2
time_period <- 200


parameters <- safir::get_parameters(
  population = pop$n,
  contact_matrix_set = squire::get_mixing_matrix(iso3c = iso3c),
  iso3c = iso3c,
  R0 = R0,
  time_period = time_period
)


dt <- 0.1



# compare to squire (the gold standard)

out <- squire::run_explicit_SEEIR_model(
  population = pop$n,
  country = "Antigua and Barbuda",
  contact_matrix_set = squire::get_mixing_matrix(iso3c = "ATG"),
  time_period = 200,
  replicates = nrep,
  day_return = TRUE,
  R0 = 2,
  dt = dt
)



# replicates of safir run
system.time(
  saf_reps <- mclapply(X = 1:nrep,FUN = function(x){
    timesteps <- parameters$time_period/dt
    variables <- safir::create_variables(pop = pop, parameters = parameters)
    events <- safir::create_events(parameters = parameters,vaccines = NULL)
    safir::attach_event_listeners(variables = variables,events = events,parameters = parameters, dt = dt)
    renderer <- individual::Render$new(timesteps)
    processes <- list(
      safir::infection_process(parameters = parameters,variables = variables,events = events,dt = dt,vaccines = NULL),
      individual::categorical_count_renderer_process(renderer, variables$state, categories = variables$states$get_categories())
    )
    safir::setup_events(parameters = parameters,events = events,variables = variables,dt = dt,vaccines = NULL)

    individual::simulation_loop(
      variables = variables,
      events = events,
      processes = processes,
      timesteps = timesteps
    )
    df <- renderer$to_dataframe()
    df$repetition <- x
    return(df)
  },mc.cores = ncores)
)

saf_reps <- do.call(rbind,saf_reps)

saf_dt <- as.data.table(saf_reps)
saf_dt[, IMild_count := IMild_count + IAsymp_count]
saf_dt[, IAsymp_count := NULL]
saf_dt <- melt(saf_dt,id.vars = c("timestep","repetition"),variable.name = "name")
saf_dt[, model := "safir"]
saf_dt[, name := gsub("(^)(\\w*)(_count)", "\\2", name)]
setnames(x = saf_dt,old = c("timestep","name","value"),new = c("t","compartment","y"))
saf_dt[, t := t * dt]
saf_dt <- saf_dt[, .(ymin = quantile(y,0.025), ymax = quantile(y,0.975), y = median(y)), by = .(t,compartment,model)]


sq_dt <- as.data.table(squire::format_output(out, unique(saf_dt$compartment)))
sq_dt[, model := "squire"]
sq_dt <- sq_dt[, .(ymin = quantile(y,0.025), ymax = quantile(y,0.975), y = median(y)), by = .(t,compartment,model)]

ggplot(data = rbind(saf_dt,sq_dt), aes(t,y,color = model)) +
  geom_line() +
  geom_ribbon(ggplot2::aes(ymin = ymin, ymax = ymax, fill = model), alpha = 0.2) +
  geom_line() +
  facet_wrap(~compartment, scales = "free")


#  recreate natl fits


# # grab the json from the data exports
# iso3c <- "GBR"
# file_path <- "https://raw.githubusercontent.com/mrc-ide/global-lmic-reports/master/"
# country <- squire::population$country[squire::population$iso3c==iso3c][1]
# json_path <- file.path(file_path,iso3c,"input_params.json")
# json <- jsonlite::read_json(json_path)
# betas <- vapply(json, "[[", FUN.VALUE = numeric(1), "beta_set")
#
# # build a 100th of UK pop for simplicity
# uk_pop <- safir:::get_population(iso3c)
# div <- 100
# uk_pop$n <- as.integer(uk_pop$n / div)
#
# # get first 365 days of fit
# uk_parameters <- safir::get_parameters(
#   population = uk_pop$n,
#   contact_matrix_set = squire::get_mixing_matrix(iso3c = iso3c),
#   iso3c = iso3c,
#   beta_set = head(betas, 365),
#   R0 = head(betas, 365),
#   dur_R = 365,
#   seeding_cases = 5,
#   tt_R0 = seq_along(head(betas, 365)),
#   time_period = length(head(betas, 365))
# )
#
# nrep <- 4
#
# system.time(
#   saf_reps_uk <- mclapply(X = 1:nrep,FUN = function(x){
#     timesteps <- uk_parameters$time_period/dt
#     variables <- safir::create_variables(pop = uk_pop, parameters = uk_parameters)
#     events <- safir::create_events(parameters = uk_parameters,vaccines = NULL)
#     safir::attach_event_listeners(variables = variables,events = events,parameters = uk_parameters, dt = dt)
#     renderer <- individual::Render$new(timesteps)
#     processes <- list(
#       safir::infection_process(parameters = uk_parameters,variables = variables,events = events,dt = dt,vaccines = NULL),
#       individual::categorical_count_renderer_process(renderer, variables$state, categories = variables$states$get_categories())
#     )
#     safir::setup_events(parameters = uk_parameters,events = events,variables = variables,dt = dt,vaccines = NULL)
#
#     individual::simulation_loop(
#       variables = variables,
#       events = events,
#       processes = processes,
#       timesteps = timesteps
#     )
#     df <- renderer$to_dataframe()
#     df$repetition <- x
#     return(df)
#   },mc.cores = ncores)
# )
#
# saf_reps_uk <- do.call(rbind,saf_reps_uk)
#
# saf <- as.data.table(saf_reps_uk)
# saf[, setdiff(names(saf), c("timestep","D_count","repetition")) := NULL]
# setnames(x = saf,old = c("timestep","D_count"),new = c("t","D"))
# saf <- melt(data = saf,id.vars = c("t","repetition"),variable.name = "compartment",value.name = "y")
# saf[, t := t*dt]
#
#
# # real deaths in the UK
# real <- data.frame(
#   "date" = as.Date(vapply(json, "[[", FUN.VALUE = character(1), "date"))[1:365],
#   "deaths" = cumsum(vapply(json, function(x){ if(!is.null(x$deaths)) x$deaths else  0}, numeric(1))[1:365]),
#   "t" = 1:365
# )
#
# ggplot(data = merge(saf, real, "t"),mapping = aes(date, y*div, group = repetition)) +
#   geom_bar(aes(date, deaths, group = "real"), fill = "red",
#                     stat = "identity", position=position_dodge(width = 0)) +
#   geom_line(color = "black") +
#   ylab("Cumulative Deaths") +
#   xlab("Date") +
#   theme_bw()
#
#
#
# # run on 4 cores
# reps <- 4
# options("mc.cores" = 2)
# saf_reps_uk <- safir::run_simulation_replicate(
#   repetitions = reps,
#   overrides = safir:::rep_list(list("pop" = uk_pop, "parameters" = uk_parameters), reps),
#   parallel = TRUE
# )
#
# # sort results and just filter to cumulative
# saf_tib <- saf_reps_uk %>%
#   dplyr::mutate(IMild_count = IMild_count + IAsymp_count) %>%
#   dplyr::select(-IAsymp_count) %>%
#   tidyr::pivot_longer(-c(timestep, repetition)) %>%
#   dplyr::mutate(name = gsub("(^)(\\w*)(_count)", "\\2", name),
#                 model = "safir") %>%
#   dplyr::rename(t = timestep, compartment = name, y = value) %>%
#   dplyr::select(t, compartment, y, model, repetition) %>%
#   dplyr::filter(compartment == "D")
#
# # real deaths in the UK
# real <- data.frame(
#   "date" = as.Date(vapply(json, "[[", FUN.VALUE = character(1), "date"))[1:365],
#   "deaths" = cumsum(vapply(json, function(x){ if(!is.null(x$deaths)) x$deaths else  0}, numeric(1))[1:365]),
#   "t" = 1:365
# )
#
# # join and plot
# dplyr::left_join(saf, real, "t") %>%
#   ggplot2::ggplot(ggplot2::aes(date, y*div, group = repetition)) +
#   ggplot2::geom_bar(ggplot2::aes(date, deaths, group = "real"), fill = "red",
#                     stat = "identity", position=ggplot2::position_dodge(width = 0)) +
#   ggplot2::geom_line(color = "black") +
#   ggplot2::ylab("Cumulative Deaths") +
#   ggplot2::xlab("Date") +
#   ggplot2::theme_bw()
