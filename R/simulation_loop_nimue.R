# --------------------------------------------------
#   Simulation loop for nimue style vaccinations
#   Sean L. Wu (slwood89@gmail.com)
#   May 2021
# --------------------------------------------------

#' @title Simulation loop for safir (nimue vaccine model)
#' @description a modification of \code{\link[individual]{simulation_loop}}
#' to deal with the presence of special values in the \code{variables} list, and
#' also adds an optional progress bar
#' @param variables a list of Variables
#' @param events a list of Events
#' @param processes a list of processes to execute on each timestep
#' @param timesteps the number of timesteps to simulate
#' @param progress dispaly a progress bar?
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
simulation_loop_nimue <- function(
  variables,
  events,
  processes,
  timesteps,
  progress = FALSE
) {
  if (timesteps <= 0) {
    stop('End timestep must be > 0')
  }
  if (progress) {
    pb <- txtProgressBar(min = 1, max = timesteps, initial = 1, style = 3)
  }
  # sim loop
  for (t in seq_len(timesteps)) {
    for (process in processes) {
      individual:::execute_any_process(process, t)
    }
    for (event in events) {
      event$.process()
    }
    # variables need special treatment, some shouldn't be updated, and age doesn't change
    variables$vaccine_states$.update()
    variables$states$.update()
    # event clocks tick
    for (event in events) {
      event$.tick()
    }
    if (progress) {
      setTxtProgressBar(pb = pb,value = t)
    }
  }
  # end loop, close out pb
  if (progress) {
    close(con = pb)
  }
}
