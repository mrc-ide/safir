library(safir)
library(data.table)
library(ggplot2)
library(parallel)

iso3c <- "ATG"
pop <- safir:::get_population(iso3c)

# use as many as you want normally.

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

var <- create_variables(pop = pop,parameters = parameters)

vaxx <- create_vaccine_variables(pop = pop)

var <- append(var,vaxx)

get_proportion_vaccinated <- function(variables, timestep, age, dose) {
  age_bset <- variables$discrete_age$get_index_of(age)
  N <- age_bset$size()
  vaccinated_bset <- variables$dose_time[[dose]]$get_index_of(a = 0,b = timestep)
  vaccinated_bset$and(age_bset)
  return( vaccinated_bset$size() / N )
}
