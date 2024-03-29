library(individual)
library(nimue)
library(squire)

tmax <- 100

iso3c <- "GBR"
pop <- safir::get_population(iso3c)
pop$n <- as.integer(pop$n / 1e3)
contact_mat <- get_mixing_matrix(iso3c = iso3c)

n <- 17 * 100
ages <- rep(1:17, each = 100)
population <- tab_bins(a = ages,nbins = 17)

# vaccine dosing
vaccine_doses <- 3
dose_period <- c(NaN, 14, 7)
vaccine_set <- rep(80,100)

# vaccine strategy
vaccine_coverage_mat <- strategy_matrix(strategy = "Elderly",max_coverage = 0.8)
next_dose_priority <- matrix(0,nrow = 2, ncol = 17)
next_dose_priority[1, 15:17] <- 1
next_dose_priority[2, 10:14] <- 1

# base parameters
parameters <- safir::get_parameters(
  population = pop$n,
  contact_matrix_set = contact_mat,
  iso3c = iso3c,
  R0 = 4,
  time_period = tmax,
  dt = 1
)

# vaccine parameters
mu_ab_list <- data.frame(name = c("Pfizer", "AstraZeneca", "Sinovac", "Moderna"),
                         mu_ab_d1 = c(13/94, 1/59, 28/164, ((185+273)/2)/321),
                         mu_ab_d2 = c(223/94, 32/59, 28/164, 654/158),
                         mu_ab_d3 = c(223/94, 32/59, 28/164, 654/158))

ab_parameters <- get_vaccine_ab_titre_parameters(vaccine = "Pfizer", max_dose = vaccine_doses,correlated = FALSE, mu_ab_list = mu_ab_list)

# combine parameters and verify
parameters <- make_vaccine_parameters(
  safir_parameters = parameters,
  vaccine_ab_parameters = ab_parameters,
  vaccine_set = vaccine_set,
  dose_period = dose_period,
  strategy_matrix = vaccine_coverage_mat,
  next_dose_priority_matrix = next_dose_priority
)

parameters$population <- population


variables <- list()
variables$discrete_age <- IntegerVariable$new(ages)

full_bset <- Bitset$new(n)$insert(1:n)
