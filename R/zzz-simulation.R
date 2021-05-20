#' @title Run the simulation
#' @description
#' The main entrypoint for the simulation. run_simulation puts together the
#' model components and runs the malaria simulation. This currently returns a
#' dataframe with the number of individuals in each state at each timestep
#'
#' @param pop population. See [squire::get_population]
#' @param parameters model parameters
#' @param dt the time step
#' @export
run_simulation_dt <- function(pop, parameters) {
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