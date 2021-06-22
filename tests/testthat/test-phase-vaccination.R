

test_that("testing get_vaccination_priority_stage", {

  parameters <- list()
  parameters$vaccine_coverage_mat <- nimue::strategy_matrix(strategy = "Elderly")
  parameters$N_age <- nrow(parameters$vaccine_coverage_mat)
  parameters$N_phase <- 3

  get_vaccination_priority_stage

})
