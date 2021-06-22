

test_that("testing get_vaccination_priority_stage", {

  parameters <- list()
  parameters$vaccine_coverage_mat <- nimue::strategy_matrix(strategy = "Elderly")
  parameters$N_age <- ncol(parameters$vaccine_coverage_mat)
  parameters$N_prioritisation_steps <- nrow(parameters$vaccine_coverage_mat)
  parameters$N_phase <- 3

  n <- 17 * 100
  ages <- rep(1:17, each = 100)

  variables <- list()
  variables$discrete_age <- IntegerVariable$new(ages)

  # prioritize 3 oldest age groups
  parameters$next_dose_priority <- matrix(0,nrow = 2, ncol = 17)
  parameters$next_dose_priority[1:2, 15:17] <- 1

  # tests for all phases and stages
  for (phase in 1:3) {
    for (stage in 0:17) {
      cat("running phase: ",phase,", stage: ",stage)

      # 0th "stage", should give stage = 1
      if (stage == 0) {

        variables$dose_time <- list()
        variables$dose_time[[1]] <- IntegerVariable$new(rep(-1,n))
        variables$dose_time[[2]] <- IntegerVariable$new(rep(-1,n))
        variables$dose_time[[3]] <- IntegerVariable$new(rep(-1,n))

        calc_stage <- get_vaccination_priority_stage(variables = variables,phase = phase,parameters = parameters)
        cat(" calc stage: ",calc_stage," --- \n")

        expect_equal(
          calc_stage,1
        )
      # final stage, should give stage = -1 (move to next)
      } else if(stage == 17) {

        variables$dose_time <- list()
        variables$dose_time[[1]] <- IntegerVariable$new(rep(-1,n))
        variables$dose_time[[2]] <- IntegerVariable$new(rep(-1,n))
        variables$dose_time[[3]] <- IntegerVariable$new(rep(-1,n))

        if (phase < 3) {
          for (p in 1:phase) {
            # groups for this dose
            variables$dose_time[[p]]$queue_update(values = 1,index = 1:n)
            variables$dose_time[[p]]$.update()
            # groups prioritized for next dose
            if (length(who_2_vaccinate_next) > 0) {
              variables$dose_time[[p + 1]]$queue_update(values = 1,index = 1:n)
              variables$dose_time[[p + 1]]$.update()
            }
          }
        } else {
          # groups for this dose
          variables$dose_time[[3]]$queue_update(values = 1,index = 1:n)
          variables$dose_time[[3]]$.update()
        }

        calc_stage <- get_vaccination_priority_stage(variables = variables,phase = phase,parameters = parameters)
        cat(" calc stage: ",calc_stage," --- \n")

        expect_equal(
          calc_stage, -1
        )
      # intermediate stages, should give stage + 1
      } else {

        variables$dose_time <- list()
        variables$dose_time[[1]] <- IntegerVariable$new(rep(-1,n))
        variables$dose_time[[2]] <- IntegerVariable$new(rep(-1,n))
        variables$dose_time[[3]] <- IntegerVariable$new(rep(-1,n))

        stage_pr_vec <- parameters$vaccine_coverage_mat[stage, ]
        who_2_vaccinate <- which(stage_pr_vec > 0)

        if (phase < 3) {
          who_2_vaccinate_next <- intersect(who_2_vaccinate, which(as.logical(parameters$next_dose_priority[phase, ])))
          for (p in 1:phase) {
            # groups for this dose
            variables$dose_time[[p]]$queue_update(values = 1,index = which(ages %in% who_2_vaccinate))
            variables$dose_time[[p]]$.update()
            # groups prioritized for next dose
            if (length(who_2_vaccinate_next) > 0) {
              variables$dose_time[[p + 1]]$queue_update(values = 1,index = which(ages %in% who_2_vaccinate_next))
              variables$dose_time[[p + 1]]$.update()
            }
          }
        } else {
          # groups for this dose
          variables$dose_time[[3]]$queue_update(values = 1,index = which(ages %in% who_2_vaccinate))
          variables$dose_time[[3]]$.update()
        }

        calc_stage <- get_vaccination_priority_stage(variables = variables,phase = phase,parameters = parameters)
        cat(" calc stage: ",calc_stage," --- \n")

        expect_equal(
          calc_stage, stage + 1
        )
      }

    }
  }

})
