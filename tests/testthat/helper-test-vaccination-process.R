library(individual)
library(nimue)

parameters <- list()
parameters$vaccine_coverage_mat <- nimue::strategy_matrix(strategy = "Elderly",max_coverage = 0.8)
parameters$N_age <- ncol(parameters$vaccine_coverage_mat)
parameters$N_prioritisation_steps <- nrow(parameters$vaccine_coverage_mat)
parameters$N_phase <- 3
parameters$dose_period <- c(NaN, 14, 7)
parameters$vaccine_set <- rep(80,100)

parameters$std10 <- 0.44
parameters$mu_ab <- c(0.14, 2.37, 1.75)
parameters$ab_50 <- 0.2 # titre relative to convalescent required to provide 50% protection from infection, on linear scale
parameters$ab_50_severe <- 0.03
parameters$k <- 2.94 # shape parameter of efficacy curve
parameters$correlated <- FALSE

parameters$next_dose_priority <- matrix(0, nrow = parameters$N_phase - 1, ncol = parameters$N_age)
parameters$next_dose_priority[1, 15:17] <- 1
parameters$next_dose_priority[2, 10:14] <- 1

n <- 17 * 100
ages <- rep(1:17, each = 100)

parameters$population <- tab_bins(a = ages,nbins = 17)

variables <- list()
variables$discrete_age <- IntegerVariable$new(ages)

full_bset <- Bitset$new(n)$insert(1:n)
