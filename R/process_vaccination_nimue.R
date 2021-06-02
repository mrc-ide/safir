vaccination_process_nimue <- function(parameters, variables, events, dt) {

  function(timestep) {

    mv <- parameters$vaccine_set[ceiling(timestep * dt)]

    if (mv > 0) {

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

      vr_den <- variables$eligible$size()
      vr <- ifelse(mv <= 0, 0, min(mv/vr_den, 1))

      variables$vr <- vr

    } else {
      variables$eligible$and(variables$empty)
      variables$vr <- 0
    }

  }
}
