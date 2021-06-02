

# this is run before anything else runs in the time step
vaccination_process_zzz <- function(parameters, variables, events) {

  function(timestep) {
    # needs to compute
    # 1. number of vaccines that can be given out today (vaccine rate * dt), because vaccine rate is in
    #    numbers of doses/day
    # 2.
  }
}

vaccination_process_nimue <- function(parameters, variables, events, dt) {

  function(timestep) {

    # calculate prioritisation step and which age groups are eligible right now
    pr <- sapply(X = 1:parameters$N_age,FUN = function(a){get_proportion_vaccinated_nimue(variables = variables,age = a)})

    vaccination_target_mat <- matrix(data = 0,nrow = parameters$N_prioritisation_steps,ncol = parameters$N_age)
    for (p in 1:parameters$N_prioritisation_steps) {
      vaccination_target_mat[p, ] <- as.integer(pr < parameters$vaccine_coverage_mat[p, ])
    }

    vaccine_target_vec <- rep(0, parameters$N_prioritisation_steps)
    for (p in 1:parameters$N_prioritisation_steps) {
      # an entire row summing to zero means that step has been completed
      vaccine_target_vec[p] <- as.integer(sum(vaccination_target_mat[p, ]) == 0)
    }
    current_index <- min(sum(vaccine_target_vec) + 1, parameters$N_prioritisation_steps)

    vaccination_target <- vaccination_target_mat[current_index, ]

    # calculate the current per-capita rate of vaccination
    unvaccinated <- variables$vaccinated$not()

    # clear out eligible
    variables$eligible$and(variables$empty)

    SER <- variables$states$get_index_of(c("S","E","R"))
    for (a in which(vaccination_target > 0)) {
      SER$and(variables$discrete_age$get_index_of(a))
    }

    # set who is eligible
    variables$eligible$or(SER)

    # needs to compute
    # 1. number of vaccines that can be given out today (vaccine rate * dt), because vaccine rate is in
    #    numbers of doses/day
    # 2.
  }
}
