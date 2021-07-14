# --------------------------------------------------
#   functions to work with scheduling/allocating vaccines
#   for the general vaccine model but with no types
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


#' @title Get prioritisation step for a specific dosing phase
#' @description Return which row of the \code{\link[nimue]{strategy_matrix}} the vaccination
#' program should be targeting for coverage. To complete a step, all \code{phase}
#' dose coverage should be > prioritisation matrix target and all \code{phase + 1}
#' dose coverage for prioritized groups should be > prioritisation matrix target.
#' If these two conditions are fulfilled for the entire \code{phase}, the function returns
#' -1 to indicate that vaccination dosing \code{phase} should be advanced.
#' @param variables a list
#' @param phase current dosing phase
#' @param parameters a list
#'
#' @export
get_vaccination_priority_stage <- function(variables, phase, parameters) {

  stopifnot(is.finite(parameters$N_phase))
  stopifnot(nrow(parameters$next_dose_priority) == parameters$N_phase - 1)
  stopifnot(ncol(parameters$next_dose_priority) == ncol(parameters$vaccine_coverage_mat))

  # not final phase
  if (phase < parameters$N_phase) {

    # pr_this_dose <- sapply(X = 1:parameters$N_age,FUN = function(a){get_proportion_vaccinated(variables = variables, age = a, dose = phase)})
    # pr_next_dose <- sapply(X = 1:parameters$N_age,FUN = function(a){get_proportion_vaccinated(variables = variables, age = a, dose = phase + 1)})

    pr_this_dose <- get_proportion_vaccinated_all_ages_cpp(variables = variables,N_age = parameters$N_age,dose = phase)
    pr_next_dose <- get_proportion_vaccinated_all_ages_cpp(variables = variables,N_age = parameters$N_age,dose = phase + 1)

    # go through prioritisation steps
    for (p in 1:parameters$N_prioritisation_steps) {
      ages_2_check <- which(parameters$vaccine_coverage_mat[p, ] > 0)
      this_dose_not_cover <- any(pr_this_dose[ages_2_check] < parameters$vaccine_coverage_mat[p, ages_2_check])

      ages_2_check_next <- which(parameters$next_dose_priority[phase, ] > 0)
      next_dose_not_cover <- any(pr_next_dose[ages_2_check_next] < parameters$vaccine_coverage_mat[p ,ages_2_check_next])

      if (this_dose_not_cover | next_dose_not_cover) {
        return(p)
      }
    }

    # if we did not return by now it means this step is complete, return -1
    return(-1)

  #  final phase: don't need to check for next dose coverage
  } else {

    # pr_this_dose <- sapply(X = 1:parameters$N_age,FUN = function(a){get_proportion_vaccinated(variables = variables, age = a, dose = phase)})

    pr_this_dose <- get_proportion_vaccinated_all_ages_cpp(variables = variables,N_age = parameters$N_age,dose = phase)

    # go through prioritisation steps
    for (p in 1:parameters$N_prioritisation_steps) {
      ages_2_check <- which(parameters$vaccine_coverage_mat[p, ] > 0)
      this_dose_not_cover <- any(pr_this_dose[ages_2_check] < parameters$vaccine_coverage_mat[p, ages_2_check])

      if (this_dose_not_cover) {
        return(p)
      }
    }

    # if we did not return by now it means this step is complete, return -1
    return(-1)

  }

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


#' @title Find individuals in each age group eligible for a dose (multi-dose, no types)
#' @description Find those individuals who have had the dose preceding \code{dose},
#' have not yet received the next one, and are beyond the \code{dose_period}.
#' This is similar to the function \code{\link[nimue]{eligable_for_second}} in the nimue package.
#' This function should only be called if simulation time is greater than the \code{dose_period}
#' @param dose which dose?
#' @param parameters a list
#' @param variables a list
#' @param t the current time step
#' @param dt size of the time step
#' @return a list containing a \code{\link[individual]{Bitset}} for each age bin
#' @export
age_group_eligible_for_dose_vaccine <- function(dose, parameters, variables, t, dt) {
  out <- replicate(n = parameters$N_age,expr = NULL,simplify = FALSE)
  if (dose > 1) {
    threshold <- t - as.integer(parameters$dose_period[dose]/dt)
    # nobody is eligible
    if (threshold < 0) {
      for (a in 1:parameters$N_age) {
        out[[a]] <- individual::Bitset$new(sum(parameters$population))
      }
      return(out)
    }
    # who has gotten the previous dose in [0,threshold]?
    had_previous_beyond_threshold <- variables$dose_time[[dose - 1]]$get_index_of(a = 0, b = threshold)
    # who has not gotten the next one?
    not_had_next <- variables$dose_time[[dose]]$get_index_of(set = -1)
    # people past the threshold and who haven't gotten the next one yet
    not_had_next$and(had_previous_beyond_threshold)
    # calculate by age
    for (a in 1:parameters$N_age) {
      bset_a <- variables$discrete_age$get_index_of(set = a)
      bset_a$and(not_had_next)
      out[[a]] <- bset_a
    }
  } else {
    not_had_first_dose <- variables$dose_time[[dose]]$get_index_of(set = -1)
    # calculate by age
    for (a in 1:parameters$N_age) {
      bset_a <- variables$discrete_age$get_index_of(set = a)
      bset_a$and(not_had_first_dose)
      out[[a]] <- bset_a
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
  # pr <- sapply(X = 1:parameters$N_age,FUN = function(a){get_proportion_vaccinated(variables = variables,age = a,dose = dose)})
  pr <- get_proportion_vaccinated_all_ages_cpp(variables = variables,N_age = parameters$N_age,dose = dose)

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


#' @title Target persons in each age group to vaccinate (multi-dose, no types)
#' @description For a given dose, find eligible persons in each age group to target for that vaccine dose.
#' If \code{vaxx_priority} is specified, it means that we are currently vaccinating one dose ahead
#' of the current vaccination phase, for those age groups prioritized to get their next dose.
#' @param dose which dose?
#' @param variables a list
#' @param parameters Model parameters
#' @param t current time step
#' @param dt size of time step
#' @param prioritisation row of the prioritisation matrix for the current step
#' @param vaxx_priority vaxx_priority matrix for age groups that should get the next dose
#' the current dose
#' @return a named list with \code{n_to_cover}: number to be vaccinated in each age group,
#' \code{eligible_age_bsets}: \code{\link[individual]{Bitset}} of eligible group in each age bin (from \code{\link{age_group_eligible_for_dose_vaccine}})
#' \code{eligible_age_counts}: size of eligible group in each age bin (from \code{\link{age_group_eligible_for_dose_vaccine}})
#' @export
target_pop <- function(dose, variables, parameters, t, dt, prioritisation, vaxx_priority = NULL) {

  # current coverage by age in this dose
  current_coverage <- sapply(X = 1:parameters$N_age,FUN = function(a){
    get_proportion_vaccinated(variables = variables,age = a,dose = dose)
  })

  # size of each age group
  age_size <- sapply(X = 1:parameters$N_age,FUN = function(a){
    variables$discrete_age$get_size_of(set = a)
  })

  # Remaining population left to cover with current dose number (dose) to reach target coverage in prioritisation step
  n_to_cover <- ceiling(pmax(0, (prioritisation - current_coverage)) * age_size)

  out <- list(
    n_to_cover = NULL,
    eligible_age_bsets = NULL,
    eligible_age_counts = NULL
  )

  # eligible group
  eligible_age_bsets <- age_group_eligible_for_dose_vaccine(dose = dose,parameters = parameters, variables = variables,t = t, dt = dt)
  eligible_age_counts <- sapply(eligible_age_bsets, function(x){x$size()})

  out$eligible_age_bsets <- eligible_age_bsets
  out$eligible_age_counts <- eligible_age_counts

  if (dose > 1) {

    if (!is.null(vaxx_priority)) {
      # giving priority vaccine doses to dose + 1
      n_to_cover <- pmin(n_to_cover, eligible_age_counts) * vaxx_priority
      out$n_to_cover <- n_to_cover
    } else {
      # normal dose; do not check vaxx_priority
      n_to_cover <- pmin(n_to_cover, eligible_age_counts)
      out$n_to_cover <- n_to_cover
    }

  } else {
    # dose 1, we can find out just from proportion with coverage the number to give
    # without checking eligibility after threshold
    out$n_to_cover <- n_to_cover
  }

  stopifnot(all(out$n_to_cover <= out$eligible_age_counts))

  return(out)
}


#' Schedule N doses to age groups based on weightings by number of people eligible to be vaccinated (multi-dose, no types)
#' @description Please make sure you are subtracting these doses from some daily total outside of this function.
#' It combines functionality of \code{\link[nimue]{assign_doses}} and \code{\link[nimue]{administer_first_dose}}/\code{\link[nimue]{administer_second_dose}}.
#' Be aware that it is possible for fewer doses to be assigned than supplied, if not enough doses are available. In general
#' \code{n_to_cover} is the ideal allocation, and it is constrained by the amount available, \code{doses}
#' @param doses Total available doses
#' @param n_to_cover Number of people eligible to be vaccinated in each age group, from \code{\link{target_pop_new}}
#' @param eligible_age_bset a list of \code{\link[individual]{Bitset}} from \code{\link{target_pop_new}}
#' @param eligible_age_counts size of eligible group in each age bin from \code{\link{target_pop_new}}
#' @param events a list
#' @param dose what dose we are administering
#' @param parameters a list
#' @param discrete_age a \code{\link[individual]{IntegerVariable}}
#' @return the number of remaining doses
#' @export
assign_doses <- function(doses, n_to_cover, eligible_age_bset, eligible_age_counts, events, dose, parameters, discrete_age) {

  stopifnot(all(n_to_cover <= eligible_age_counts))

  leftover_doses <- doses

  # check if any vaccines are scheduled at all
  if (any(n_to_cover > 1 )) {

    # no dose scarcity
    if (sum(n_to_cover) < leftover_doses) {

      # loop over ages
      for (a in 1:parameters$N_age) {

        # allocated doses < eligible persons, schedule subset
        if (eligible_age_counts[a] != n_to_cover[a]) {

          num_to_retain <- n_to_cover[a]
          n <- eligible_age_counts[a]
          to_keep <- sample.int(n = n,size = num_to_retain,replace = FALSE)

          to_schedule <- filter_bitset(eligible_age_bset[[a]], to_keep)
          to_schedule$and(events$scheduled_dose[[dose]]$get_scheduled()$not())
          events$scheduled_dose[[dose]]$schedule(target = to_schedule, delay = 0)

          leftover_doses <- leftover_doses - to_schedule$size()

        # schedule everyone
        } else {
          to_schedule <- eligible_age_bset[[a]]
          to_schedule$and(events$scheduled_dose[[dose]]$get_scheduled()$not())
          events$scheduled_dose[[dose]]$schedule(target = to_schedule, delay = 0)

          leftover_doses <- leftover_doses - to_schedule$size()
        }

      } # end age loop

    # dose scarcity; need to allocate proportional to group size
    } else {

      # people already getting a dose, by age
      already_scheduled <- tab_bins(a = discrete_age$get_values(events$scheduled_dose[[dose]]$get_scheduled()), nbins = length(n_to_cover))
      # don't want to assign vaccines to those already scheduled for one
      n_to_cover <- pmax(0, n_to_cover - already_scheduled)
      group_weights <- n_to_cover / sum(n_to_cover)
      # make sure assigned vaccines never exceed the number of leftover doses
      assigned <- floor(min(sum(n_to_cover), leftover_doses) * group_weights)
      if(sum(assigned) != leftover_doses){
        assigned <- assigned + (rank(-group_weights, ties.method = "last") <= (leftover_doses %% length(group_weights)))
      }

      # loop over ages
      for (a in 1:parameters$N_age) {

        num_to_retain <- assigned[a]
        n <- eligible_age_counts[a]
        to_keep <- sample.int(n = n,size = num_to_retain,replace = FALSE)

        to_schedule <- filter_bitset(eligible_age_bset[[a]], to_keep)
        to_schedule$and(events$scheduled_dose[[dose]]$get_scheduled()$not())
        events$scheduled_dose[[dose]]$schedule(target = to_schedule, delay = 0)

        leftover_doses <- leftover_doses - to_schedule$size()

      } # end age loop

    } # end dose scarcity

  } # end allocation

  stopifnot(leftover_doses >= 0) # assertion
  return(leftover_doses)
}

