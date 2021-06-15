# --------------------------------------------------
#   functions to work with scheduling/allocating vaccines
#   for the general vaccine model
#   Sean L. Wu (slwood89@gmail.com)
#   June 2021
# --------------------------------------------------

#' @title Get proportion of an age group that is vaccinated
#' @description Get proportion of an age group that has received a particular vaccine dose
#' by this timestep. This is similar to the function \code{\link[nimue]{coverage}} in the nimue package.
#' @param variables a list
#' @param age an age group (can be a vector of Multiple age groups)
#' @param dose which dose to get proportion of age group who has received it
#'
#' @export
get_proportion_vaccinated <- function(variables, age, dose) {
  age_bset <- variables$discrete_age$get_index_of(age)
  N <- age_bset$size()
  vaccinated_bset <- variables$dose_type[[dose]]$get_index_of("UNVACC") # haven't gotten this dose
  vaccinated_bset <- vaccinated_bset$not() # complement = vaccinated people
  vaccinated_bset$and(age_bset)
  return( vaccinated_bset$size() / N )
}


#' @title Get proportion of an age group that is vaccinated, by type
#' @description Get proportion of an age group that has received a particular vaccine type and dose
#' by this timestep. This is similar to the function \code{\link[nimue]{coverage}} in the nimue package.
#' @param variables a list
#' @param age an age group (can be a vector of Multiple age groups)
#' @param type type of vaccine
#' @param dose which dose to get proportion of age group who has received it
#'
#' @export
get_proportion_vaccinated_type <- function(variables, age, type, dose) {
  # stopifnot(type %in% variables$dose_type[[dose]]$get_categories())
  age_bset <- variables$discrete_age$get_index_of(age)
  N <- age_bset$size()
  type_bset <- variables$dose_type[[dose]]$get_index_of(type)
  vaccinated_bset <- variables$dose_type[[dose]]$get_index_of("UNVACC") # haven't gotten this dose
  vaccinated_bset <- vaccinated_bset$not() # complement = vaccinated people
  vaccinated_bset$and(age_bset)$and(type_bset)
  return( vaccinated_bset$size() / N )
}


#' @title Identity those persons eligible for a dose
#' @description Find those individuals who have had the dose preceding \code{dose},
#' have not yet received the next one, and are beyond the \code{dose_period}.
#' This is similar to the function \code{\link[nimue]{eligable_for_second}} in the nimue package.
#' This function should only be called if simulation time is greater than the \code{dose_period}
#' @param dose which dose?
#' @param dose_period days between \code{dose} and the previous dose (please make sure \code{dose / dt} produces an integer number of timesteps)
#' @param variables a list
#' @param t the current time step
#' @param dt size of the time step
#' @return an \code{\link[individual]{Bitset}}
#' @export
eligible_for_dose_vaccine <- function(dose, dose_period, variables, t, dt) {
  if (dose > 1) {
    # who has gotten the previous dose? (with correction for dt < 1)
    had_previous_beyond_threshold <- variables$dose_time[[dose - 1]]$get_index_of(a = 0, b = t - as.integer(dose_period/dt))
    # who has not gotten the next one?
    not_had_next <- variables$dose_time[[dose]]$get_index_of(set = -1)
    # return people past the threshold and who haven't gotten the next one yet
    return(not_had_next$and(had_previous_beyond_threshold))
  } else {
    not_had_first_dose <- variables$dose_time[[dose]]$get_index_of(set = -1)
    return(not_had_first_dose)
  }
}


#' @title Get current prioritization step for a specific dose
#' @param variables a list
#' @param parameters Model parameters
#' @param dose which dose?
#' @export
get_current_prioritization_step <- function(variables, parameters, dose) {
  # calculate prioritisation step and which age groups are eligible right now
  pr <- sapply(X = 1:parameters$N_age,FUN = function(a){get_proportion_vaccinated(variables = variables,age = a,dose = dose)})

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

  return(current_index)
}




# eligible_for_dose_vaccine_type same as above but with types

target_assign_pop_vaccine <- function(t, dt, doses, phase, variables, parameters, events) {

  # what step of prioritization matrix are we on?
  step <- get_current_prioritization_step(..., dose = phase)
  prioritisation <- parameters$vaccine_coverage_mat[step, ]

  # current coverage by age in this phase
  current_coverage <- sapply(X = 1:parameters$N_age,FUN = function(a){
    get_proportion_vaccinated(variables = variables,age = a,dose = phase)
  })

  # size of each age group
  age_size <- sapply(X = 1:parameters$N_age,FUN = function(a){
    variables$discrete_age$get_size_of(set = a)
  })

  # how many people need to get vaccinated in each age group to hit coverage
  age_vaxx <- ceiling(pmax(0, (prioritisation - current_coverage)) * age_size)

  # distribute vaccines

  # everyone who is eligible
  eligible <- eligible_for_dose_vaccine(dose = phase,dose_period = parameters$dose_period[dose],variables = variables,t = t, dt = dt)

  eligible_age_counts <- rep(0, parameters$N_age)
  eligible_age_bset <- replicate(n = parameters$N_age,expr = NULL,simplify = FALSE)
  for (a in 1:parameters$N_age) {
    # who is eligible and in this age group?
    eligible_age_bset[a] <- variables$discrete_age$get_index_of(set = a)
    eligible_age_bset[a]$and(eligible)
    # put size of group in a vector
    eligible_age_counts[a] <- eligible_age_bset[a]$size()
  }

  # everybody gets a dose
  if (sum(eligible_age_counts) <= doses) {
    # queue a vaccine for everyone (queue event for everyone)
    # now use the vector assigned to send out vaccines
    for (a in 1:parameters$N_age) {
      # event$thing$schedule(eligible_age_bset[a])
    }
  } else {
    # dose scarcity, need to distribute doses with weights proportional to group size
    group_weights <- eligible_age_counts / sum(eligible_age_counts)
    assigned <- floor(doses * group_weights)
    if(sum(assigned) != doses){
      assigned <- assigned + (rank(-group_weights, ties.method = "last") <= (doses %% length(group_weights)))
    }
    # now use the vector assigned to send out vaccines
    for (a in 1:parameters$N_age) {
      # event$thing$schedule(eligible_age_bset[a]$sample_indices(num_to_retain = assigned[a]))
    }
  }

}
