# --------------------------------------------------
#   functions to work with scheduling/allocating vaccines
#   for the general vaccine model
#   Sean L. Wu (slwood89@gmail.com)
#   June 2021
# --------------------------------------------------


#' @title Get proportion of an age group that is vaccinated (multi-dose, no types)
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
  vaccinated_bset <- variables$dose_time[[dose]]$get_index_of(set = -1) # haven't gotten this dose
  vaccinated_bset <- vaccinated_bset$not() # complement = vaccinated people
  vaccinated_bset$and(age_bset)
  return( vaccinated_bset$size() / N )
}


#' @title Get proportion of an age group that is vaccinated (multi-dose, with types)
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
  vaccinated_bset <- variables$dose_time[[dose]]$get_index_of(set = -1) # haven't gotten this dose
  vaccinated_bset <- vaccinated_bset$not() # complement = vaccinated people
  vaccinated_bset$and(age_bset)$and(type_bset)
  return( vaccinated_bset$size() / N )
}


#' @title Identity those persons eligible for a dose (multi-dose, no types)
#' @description Find those individuals who have had the dose preceding \code{dose},
#' have not yet received the next one, and are beyond the \code{dose_period}.
#' This is similar to the function \code{\link[nimue]{eligable_for_second}} in the nimue package.
#' This function should only be called if simulation time is greater than the \code{dose_period}
#' @param dose which dose?
#' @param parameters a list
#' @param variables a list
#' @param t the current time step
#' @param dt size of the time step
#' @return an \code{\link[individual]{Bitset}}
#' @export
eligible_for_dose_vaccine <- function(dose, parameters, variables, t, dt) {
  if (dose > 1) {
    # who has gotten the previous dose? (with correction for dt < 1)
    threshold <- t - as.integer(parameters$dose_period[dose]/dt)
    if (threshold < 0) {
      return(individual::Bitset$new(sum(parameters$population)))
    }
    had_previous_beyond_threshold <- variables$dose_time[[dose - 1]]$get_index_of(a = 0, b = threshold)
    # who has not gotten the next one?
    not_had_next <- variables$dose_time[[dose]]$get_index_of(set = -1)
    # people past the threshold and who haven't gotten the next one yet
    return(not_had_next$and(had_previous_beyond_threshold))
  } else {
    not_had_first_dose <- variables$dose_time[[dose]]$get_index_of(set = -1)
    return(not_had_first_dose)
  }
}

#' @title Find number of individuals in each age group eligible for a dose (multi-dose, no types)
#' @description Find those individuals who have had the dose preceding \code{dose},
#' have not yet received the next one, and are beyond the \code{dose_period}.
#' This is similar to the function \code{\link[nimue]{eligable_for_second}} in the nimue package.
#' This function should only be called if simulation time is greater than the \code{dose_period}
#' @param dose which dose?
#' @param parameters a list
#' @param variables a list
#' @param t the current time step
#' @param dt size of the time step
#' @return a vector of population sizes
#' @export
age_group_eligible_for_dose_vaccine <- function(dose, parameters, variables, t, dt) {
  out <- rep(0, parameters$N_age)
  if (dose > 1) {
    # who has gotten the previous dose? (with correction for dt < 1)
    threshold <- t - as.integer(parameters$dose_period[dose]/dt)
    if (threshold < 0) {
      return(out)
    }
    had_previous_beyond_threshold <- variables$dose_time[[dose - 1]]$get_index_of(a = 0, b = threshold)
    # who has not gotten the next one?
    not_had_next <- variables$dose_time[[dose]]$get_index_of(set = -1)
    # people past the threshold and who haven't gotten the next one yet
    not_had_next$and(had_previous_beyond_threshold)
    # calculate by age
    for (a in 1:parameters$N_age) {
      bset_a <- variables$discrete_age$get_index_of(set = a)
      bset_a$and(not_had_next)
      out[a] <- bset_a$size()
    }
  } else {
    not_had_first_dose <- variables$dose_time[[dose]]$get_index_of(set = -1)
    # calculate by age
    for (a in 1:parameters$N_age) {
      bset_a <- variables$discrete_age$get_index_of(set = a)
      bset_a$and(not_had_first_dose)
      out[a] <- bset_a$size()
    }
  }
  return(out)
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


#' @title Target number of individual in each age group to vaccinate (multi-dose, no types)
#' @param phase which dose?
#' @param variables a list
#' @param parameters Model parameters
#' @param t current time step
#' @param dt size of time step
#' @param prioritisation row of the prioritisation matrix for the current step
#' @param vaxx_priority row of the vaxx_priority matrix for age groups that should get the next dose phase while still on
#' the current phase
#' @export
target_pop <- function(phase, variables, parameters, t, dt, prioritisation, vaxx_priority = NULL) {

  if (phase == variables$phase & phase == parameters$N_phase & !is.null(vaxx_priority)) {
    stop("on final phase no vaxx_priority for next phase should NULL")
  }
  if (variables$phase == phase & !is.null(vaxx_priority)) {
    stop("if giving vaccines to the current phase (not prioritized for phase + 1), vaxx_priority should be NULL")
  }

  # current coverage by age in this phase
  current_coverage <- sapply(X = 1:parameters$N_age,FUN = function(a){
    get_proportion_vaccinated(variables = variables,age = a,dose = phase)
  })

  # size of each age group
  age_size <- sapply(X = 1:parameters$N_age,FUN = function(a){
    variables$discrete_age$get_size_of(set = a)
  })

  # Remaining population left to cover with current dose number (phase) to reach target coverage in prioritisation step
  n_to_cover <- ceiling(pmax(0, (prioritisation - current_coverage)) * age_size)

  if (phase > 1) {

    # eligible group
    eligible <- age_group_eligible_for_dose_vaccine(dose = phase,parameters = parameters, variables = variables,t = t, dt = dt)

    if (!is.null(vaxx_priority)) {
      # giving priority vaccine doses to phase + 1
      n_to_cover <- pmin(n_to_cover, eligible) * vaxx_priority
      return(n_to_cover)
    } else {
      # normal phase; do not check vaxx_priority
      n_to_cover <- pmin(n_to_cover, eligible)
      return(n_to_cover)
    }

  } else {
    # phase 1, we can find out just from proportion with coverage the number to give
    return(n_to_cover)
  }

}


#' Assign N doses to age groups based on weightings by number of people eligible to be vaccinated (multi-dose, no types)
#' @description Please make sure you are subtracting these doses from some daily total outside of this function.
#' It combines functionality of \code{\link[nimue]{assign_doses}} and \code{\link[nimue]{administer_first_dose}}/\code{\link[nimue]{administer_second_dose}}.
#' @param t current time step
#' @param dt size of time step
#' @param doses Total available doses
#' @param n_to_cover Number of people eligible to be vaccinated in each age group, from \code{\link{target_pop}}
#' @param variables a list
#' @param events a list
#' @param phase vaccination phase (which doses are we administering)
#' @param parameters a list
#'
#' @export
assign_doses <- function(t, dt, doses, n_to_cover, variables, events, phase, parameters) {

  # eligible people by age
  eligible <- eligible_for_dose_vaccine(dose = phase,parameters = parameters,variables = variables,t = t, dt = dt)
  eligible_age_counts <- rep(0, parameters$N_age)
  eligible_age_bset <- replicate(n = parameters$N_age,expr = NULL,simplify = FALSE)
  for (a in 1:parameters$N_age) {
    # who is eligible and in this age group?
    eligible_age_bset[[a]] <- variables$discrete_age$get_index_of(set = a)
    eligible_age_bset[[a]]$and(eligible)
    # put size of group in a vector
    eligible_age_counts[a] <- eligible_age_bset[[a]]$size()
  }

  # no dose scarcity
  if (sum(n_to_cover) <= doses) {

    # queue a vaccine for everyone (queue event for everyone)
    # now use the vector assigned to send out vaccines
    for (a in 1:parameters$N_age) {
      stopifnot(eligible_age_counts[a] != n_to_cover[a]) # take this out when done debugging
      # event$thing$schedule(eligible_age_bset[a]) SCHEDULE IT
      events$scheduled_dose[[phase]]$schedule(target = eligible_age_bset[[a]], delay = 0)
    }

  } else {
    # dose scarcity; need to allocate proportional to group size
    group_weights <- eligible_age_counts / sum(eligible_age_counts)
    assigned <- floor(doses * group_weights)
    if(sum(assigned) != doses){
      assigned <- assigned + (rank(-group_weights, ties.method = "last") <= (doses %% length(group_weights)))
    }
    # now use the vector assigned to send out vaccines
    for (a in 1:parameters$N_age) {

      num_to_retain <- assigned[a]
      n <- eligible_age_counts[a]
      to_keep <- sample.int(n = n,size = num_to_retain,replace = FALSE)

      events$scheduled_dose[[phase]]$schedule(target = filter_bitset(eligible_age_bset[[a]], to_keep), delay = 0)

      # eligible_age_bset[[a]]$remove(to_remove)
      # events$scheduled_dose[[phase]]$schedule(target = eligible_age_bset[[a]], delay = 0)
    }

  }

}


# # eligible_for_dose_vaccine_type same as above but with types
#
# target_assign_pop_vaccine <- function(t, dt, doses, phase, variables, events, parameters, events) {
#
#   # what step of prioritization matrix are we on?
#   step <- get_current_prioritization_step(variables = variables,parameters = parameters, dose = phase)
#   prioritisation <- parameters$vaccine_coverage_mat[step, ]
#
#   # current coverage by age in this phase
#   current_coverage <- sapply(X = 1:parameters$N_age,FUN = function(a){
#     get_proportion_vaccinated(variables = variables,age = a,dose = phase)
#   })
#
#   # size of each age group
#   age_size <- sapply(X = 1:parameters$N_age,FUN = function(a){
#     variables$discrete_age$get_size_of(set = a)
#   })
#
#   # Remaining population left to cover with current dose number (phase) to reach target coverage in prioritisation step
#   n_to_cover <- ceiling(pmax(0, (prioritisation - current_coverage)) * age_size)
#
#   # 1. distribute doses of this phase to age groups who are eligible and in this step
#   eligible <- eligible_for_dose_vaccine(dose = phase,dose_period = parameters$dose_period[dose],variables = variables,t = t, dt = dt)
#
#
#   # 2. give out remaining doses to phase + 1 age groups who are eligible, in this step, and prioritized for phase + 1
#
#   # 3. give out remaining doses by going through rest of the prioritization matrix
#
#
#   # # everyone who is eligible
#   # eligible <- eligible_for_dose_vaccine(dose = phase,dose_period = parameters$dose_period[dose],variables = variables,t = t, dt = dt)
#   #
#   # eligible_age_counts <- rep(0, parameters$N_age)
#   # eligible_age_bset <- replicate(n = parameters$N_age,expr = NULL,simplify = FALSE)
#   # for (a in 1:parameters$N_age) {
#   #   # who is eligible and in this age group?
#   #   eligible_age_bset[a] <- variables$discrete_age$get_index_of(set = a)
#   #   eligible_age_bset[a]$and(eligible)
#   #   # put size of group in a vector
#   #   eligible_age_counts[a] <- eligible_age_bset[a]$size()
#   # }
#   #
#   # # everybody gets a dose
#   # if (sum(eligible_age_counts) <= doses) {
#   #   # queue a vaccine for everyone (queue event for everyone)
#   #   # now use the vector assigned to send out vaccines
#   #   for (a in 1:parameters$N_age) {
#   #     # event$thing$schedule(eligible_age_bset[a])
#   #   }
#   # } else {
#   #   # dose scarcity, need to distribute doses with weights proportional to group size
#   #   group_weights <- eligible_age_counts / sum(eligible_age_counts)
#   #   assigned <- floor(doses * group_weights)
#   #   if(sum(assigned) != doses){
#   #     assigned <- assigned + (rank(-group_weights, ties.method = "last") <= (doses %% length(group_weights)))
#   #   }
#   #   # now use the vector assigned to send out vaccines
#   #   for (a in 1:parameters$N_age) {
#   #     # event$thing$schedule(eligible_age_bset[a]$sample_indices(num_to_retain = assigned[a]))
#   #   }
#   # }
#
# }


