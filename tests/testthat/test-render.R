test_that("cross tab for compartments/age works", {

  SIR <- c("S", "I", "R")
  parameters <- list(N_age = 6)

  ages <- sample.int(n = parameters$N_age, size = 100, replace = TRUE)
  compartments <- sample(x = SIR, size = 100, replace = TRUE)
  sir_to_int <- c("S" = 1, "I" = 2, "R" = 3)

  comp1 <- matrix(data = 0, nrow = parameters$N_age, ncol = 3)
  for (i in 1:100) {
    comp1[ages[i], sir_to_int[[compartments[i]]]] <- comp1[ages[i], sir_to_int[[compartments[i]]]] + 1
  }

  expect_equal(colSums(comp1), as.vector(table(compartments)[SIR]))
  expect_equal(rowSums(comp1), as.vector(table(ages)[as.character(1:parameters$N_age)]))

  comp_variable <- CategoricalVariable$new(categories = SIR, initial_values = compartments)
  age_variable <- IntegerVariable$new(initial_values = ages)

  comp_render <- Render$new(timesteps = 1)

  render_proc <- compartments_age_render_process_daily(renderer = comp_render, age = age_variable, compartments = comp_variable, parameters = parameters, dt = 1)
  render_proc(t = 0.25) # no effect unless it hits one day
  render_out <- comp_render$to_dataframe()

  expect_equal(ncol(render_out), 1)

  render_proc(t = 1)
  render_out <- comp_render$to_dataframe()
  render_out <- render_out[, -1]

  comp2 <- matrix(data = 0, nrow = parameters$N_age, ncol = 3)

  for (i in 1:ncol(render_out)) {
    c_i <- strsplit(x = colnames(render_out)[i], split = "_")[[1]][2]
    a_i <- as.integer(strsplit(x = colnames(render_out)[i], split = "_")[[1]][4])
    comp2[a_i , sir_to_int[[c_i]]] <- render_out[1, i]
  }

  expect_equal(comp1, comp2)
})
