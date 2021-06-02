

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

    vaccination_target_mat <- matrix(data = 0,nrow = parameters$N_prioritisation_steps)


    # needs to compute
    # 1. number of vaccines that can be given out today (vaccine rate * dt), because vaccine rate is in
    #    numbers of doses/day
    # 2.
  }
}
