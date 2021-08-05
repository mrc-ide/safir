
# if this test is failing, uncomment the "cat" lines to find out which stage/phase
# is causing failure, then debug from there
test_that("testing get_vaccination_priority_stage for proper prioritization matrix stage and vaccination phase", {

  n <- 17 * 100
  ages <- rep(1:17, each = 100)

  parameters <- list()
  parameters$population <- n
  parameters$vaccine_coverage_mat <- nimue::strategy_matrix(strategy = "Elderly",max_coverage = 0.8)
  parameters$N_age <- ncol(parameters$vaccine_coverage_mat)
  parameters$N_prioritisation_steps <- nrow(parameters$vaccine_coverage_mat)
  parameters$N_phase <- 3
  parameters$population <- tab_bins(a = ages,nbins = 17)
  parameters$std10 <- 0.44
  parameters$mu_ab <- c(0.14, 2.37, 1.75)
  parameters$ab_50 <- 0.2 # titre relative to convalescent required to provide 50% protection from infection, on linear scale
  parameters$ab_50_severe <- 0.03
  parameters$k <- 2.94 # shape parameter of efficacy curve
  parameters$correlated <- FALSE

  events <- list(
    scheduled_dose = replicate(n = parameters$N_phase,expr = {TargetedEvent$new(population_size = n)})
  )

  variables <- list()
  variables$discrete_age <- IntegerVariable$new(ages)

  full_bset <- Bitset$new(n)$insert(1:n)

  # prioritize 3 oldest age groups
  parameters$next_dose_priority <- matrix(0,nrow = 2, ncol = 17)
  parameters$next_dose_priority[1:2, 15:17] <- 1

  # tests for all phases and stages
  for (phase in 1:3) {
    for (stage in 0:17) {
      # cat("running phase: ",phase,", stage: ",stage)

      # 0th "stage", should give stage = 1
      if (stage == 0) {

        variables <- create_vaccine_variables(variables = variables,parameters = parameters)

        calc_stage <- get_vaccination_priority_stage(variables = variables,events = events,phase = phase,parameters = parameters)
        # cat(" calc stage: ",calc_stage," --- \n")

        expect_equal(
          calc_stage,1
        )
      # final stage, should give stage = -1 (move to next)
      } else if(stage == 17) {

        variables <- create_vaccine_variables(variables = variables,parameters = parameters)

        if (phase < 3) {
          for (p in 1:phase) {
            # groups for this dose
            schedule_dose_vaccine(timestep = 1,variables = variables,target = full_bset,dose = p, parameters = parameters)
            # groups prioritized for next dose
            if (length(who_2_vaccinate_next) > 0) {
              schedule_dose_vaccine(timestep = 1,variables = variables,target = full_bset,dose = p + 1, parameters = parameters)
            }
            update_vaccine_variables(variables = variables)
          }
        } else {
          # groups for this dose
          schedule_dose_vaccine(timestep = 1,variables = variables,target = full_bset,dose = 3, parameters = parameters)
          update_vaccine_variables(variables = variables)
        }

        calc_stage <- get_vaccination_priority_stage(variables = variables,events = events,phase = phase,parameters = parameters)
        # cat(" calc stage: ",calc_stage," --- \n")

        expect_equal(
          calc_stage, -1
        )
      # intermediate stages, should give stage + 1
      } else {

        variables <- create_vaccine_variables(variables = variables,parameters = parameters)

        stage_pr_vec <- parameters$vaccine_coverage_mat[stage, ]
        who_2_vaccinate <- which(stage_pr_vec > 0)

        if (phase < 3) {
          who_2_vaccinate_next <- intersect(who_2_vaccinate, which(as.logical(parameters$next_dose_priority[phase, ])))
          for (p in 1:phase) {
            # groups for this dose
            schedule_dose_vaccine(timestep = 1,variables = variables,target = filter_bitset(full_bset, which(ages %in% who_2_vaccinate)),dose = p, parameters = parameters)

            # groups prioritized for next dose
            if (length(who_2_vaccinate_next) > 0) {
              schedule_dose_vaccine(timestep = 1,variables = variables,target = filter_bitset(full_bset, which(ages %in% who_2_vaccinate_next)),dose = p + 1, parameters = parameters)
            }

            update_vaccine_variables(variables = variables)
          }
        } else {
          # groups for this dose
          schedule_dose_vaccine(timestep = 1,variables = variables,target = filter_bitset(full_bset, which(ages %in% who_2_vaccinate)),dose = 3, parameters = parameters)
          update_vaccine_variables(variables = variables)

        }

        calc_stage <- get_vaccination_priority_stage(variables = variables,events = events,phase = phase,parameters = parameters)
        # cat(" calc stage: ",calc_stage," --- \n")

        expect_equal(
          calc_stage, stage + 1
        )
      }

    }
  }

})




# test_that("testing get_vaccination_priority_stage for proper prioritization matrix stage and vaccination phase, randomized assignment of vaccines", {
#
#   n <- 17 * 100
#   ages <- rep(1:17, each = 100)
#
#   parameters <- list()
#   parameters$vaccine_coverage_mat <- nimue::strategy_matrix(strategy = "Elderly",max_coverage = 0.8)
#   parameters$N_age <- ncol(parameters$vaccine_coverage_mat)
#   parameters$N_prioritisation_steps <- nrow(parameters$vaccine_coverage_mat)
#   parameters$N_phase <- 3
#   parameters$population <- tab_bins(a = ages,nbins = 17)
#
#   events <- list(
#     scheduled_dose = replicate(n = parameters$N_phase,expr = {TargetedEvent$new(population_size = n)})
#   )
#
#   variables <- list()
#   variables$discrete_age <- IntegerVariable$new(ages)
#
#   full_bset <- Bitset$new(n)$insert(1:n)
#
#   # prioritize 3 oldest age groups
#   parameters$next_dose_priority <- matrix(0,nrow = 2, ncol = 17)
#   parameters$next_dose_priority[1:2, 15:17] <- 1
#
#   # tests for all phases and stages
#   for (phase in 1:3) {
#     for (stage in 0:17) {
#       cat("running phase: ",phase,", stage: ",stage)
#
#       # 0th "stage", should give stage = 1
#       if (stage == 0) {
#
#         variables <- create_vaccine_variables(variables = variables,pop = n,max_dose = parameters$N_phase)
#
#         calc_stage <- get_vaccination_priority_stage(variables = variables,events = events,phase = phase,parameters = parameters)
#         cat(" calc stage: ",calc_stage," --- \n")
#
#         expect_equal(
#           calc_stage,1
#         )
#         # final stage, should give stage = -1 (move to next)
#       } else if(stage == 17) {
#
#         variables <- create_vaccine_variables(variables = variables,pop = n,max_dose = parameters$N_phase)
#
#         if (phase < 3) {
#           for (p in 1:phase) {
#             # groups for this dose
#             schedule_dose_vaccine(timestep = 1,variables = variables,target = full_bset,dose = p)
#             # groups prioritized for next dose
#             if (length(who_2_vaccinate_next) > 0) {
#               schedule_dose_vaccine(timestep = 1,variables = variables,target = full_bset,dose = p + 1)
#             }
#             update_vaccine_variables(variables = variables)
#           }
#         } else {
#           # groups for this dose
#           schedule_dose_vaccine(timestep = 1,variables = variables,target = full_bset,dose = 3)
#           update_vaccine_variables(variables = variables)
#         }
#
#         calc_stage <- get_vaccination_priority_stage(variables = variables,events = events,phase = phase,parameters = parameters)
#         cat(" calc stage: ",calc_stage," --- \n")
#
#         expect_equal(
#           calc_stage, -1
#         )
#         # intermediate stages, should give stage + 1
#       } else {
#
#         variables <- create_vaccine_variables(variables = variables,pop = n,max_dose = parameters$N_phase)
#
#         stage_pr_vec <- parameters$vaccine_coverage_mat[stage, ]
#         who_2_vaccinate <- which(stage_pr_vec > 0)
#
#         if (phase < 3) {
#           who_2_vaccinate_next <- intersect(who_2_vaccinate, which(as.logical(parameters$next_dose_priority[phase, ])))
#           for (p in 1:phase) {
#             vax_group <- which(ages %in% who_2_vaccinate)
#             vax_group <- vax_group[sample.int(n = length(vax_group),size = ceiling(length(vax_group)*0.9),replace = FALSE)]
#             # groups for this dose
#             schedule_dose_vaccine(timestep = 1,variables = variables,target = filter_bitset(full_bset, vax_group),dose = p)
#             # groups prioritized for next dose
#             if (length(who_2_vaccinate_next) > 0) {
#               vax_group <- which(ages %in% who_2_vaccinate_next)
#               vax_group <- vax_group[sample.int(n = length(vax_group),size = ceiling(length(vax_group)*0.9),replace = FALSE)]
#               schedule_dose_vaccine(timestep = 1,variables = variables,target = filter_bitset(full_bset, vax_group),dose = p + 1)
#             }
#             update_vaccine_variables(variables = variables)
#           }
#         } else {
#           vax_group <- which(ages %in% who_2_vaccinate)
#           vax_group <- vax_group[sample.int(n = length(vax_group),size = ceiling(length(vax_group)*0.9),replace = FALSE)]
#
#           # groups for this dose
#           schedule_dose_vaccine(timestep = 1,variables = variables,target = filter_bitset(full_bset, vax_group),dose = 3)
#           update_vaccine_variables(variables = variables)
#         }
#
#         calc_stage <- get_vaccination_priority_stage(variables = variables,events = events,phase = phase,parameters = parameters)
#         cat(" calc stage: ",calc_stage," --- \n")
#
#         expect_equal(
#           calc_stage, stage + 1
#         )
#       }
#
#     }
#   }
#
# })
