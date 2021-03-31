#' @title Run the simulation
#' @description
#' The main entrypoint for the simulation. run_simulation puts together the
#' model components and runs the malaria simulation. This currently returns a
#' dataframe with the number of individuals in each state at each timestep
#'
#' @param pop population. See [squire::get_population]
#' @param parameters model parameters
#' @export
run_simulation <- function(pop, parameters) {
  timesteps <- parameters$time_period
  variables <- create_variables(pop, parameters)
  events <- create_events(parameters)
  attach_event_listeners(variables, events, parameters)
  renderer <- individual::Render$new(timesteps)
  processes = create_processes(
    events,
    variables,
    parameters,
    renderer
  )
  create_setup_process(
    parameters,
    events,
    variables
  )

  individual::simulation_loop(variables = variables,
                              events = events,
                              processes = processes,
                              timesteps = timesteps)
  output <- renderer$to_dataframe()
  return(output)
}

#' @title Run the simulation with repetitions
#'
#' @param repetitions n times to run the simulation
#' @param parallel Should it be run in parallel
#' @param overrides a named list of parameters to use instead of defaults
#' @return data frame for runs
#' @export
run_simulation_replicate <- function(
  repetitions,
  overrides = list(),
  parallel = FALSE
) {

  # check the overrides
  assert_in(c("pop", "parameters"), names(overrides))
  assert_logical(parallel)

  # check lengths
  assert_eq(repetitions, length(overrides$pop))
  assert_eq(repetitions, length(overrides$parameters))

  if (parallel) {
    fapply <- parallel::mclapply
  } else {
    fapply <- lapply
  }

  # Run mutiple sims running sequentially only
  dfs <- fapply(

    seq(repetitions),
    function(repetition) {
      df <- run_simulation(pop = overrides$pop[[repetition]],
                           parameters = overrides$parameters[[repetition]])
      df$repetition <- repetition
      df
    }

  )

  do.call("rbind", dfs)

}
