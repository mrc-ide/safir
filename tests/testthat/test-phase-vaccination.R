

# test_that("testing get_vaccination_priority_stage", {
#
#   parameters <- list()
#   parameters$vaccine_coverage_mat <- nimue::strategy_matrix(strategy = "Elderly")
#   parameters$N_age <- nrow(parameters$vaccine_coverage_mat)
#   parameters$N_phase <- 3
#
#   # prioritize 3 oldest age groups
#   parameters$next_dose_priority <- matrix(0,nrow = 2, ncol = 17)
#   parameters$next_dose_priority[1:2, 15:17] <- 1
#
#   get_vaccination_priority_stage
#
# })
