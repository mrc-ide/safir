# --------------------------------------------------
#   Simulation loop for variants
#   Sean L. Wu (slwood89@gmail.com)
#   August 2021
# --------------------------------------------------

#' @title Simulation loop for safir (variants of concern, multi-dose, no types vaccine model)
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
simulation_loop_voc <- function(
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

  stopifnot("voc" %in% names(variables))

  events_2_loop <- which(names(events) != "scheduled_dose")

  # sim loop
  for (t in seq_len(timesteps)) {

    # call processes
    for (process in processes) {
      individual:::execute_any_process(process, t)
    }

    # call events
    for (event in events[events_2_loop]) {
      event$.process()
    }
    for (event in events$scheduled_dose) {
      event$.process()
    }

    # state updates
    update_vaccine_variables(variables = variables)
    update_voc_variables(variables = variables)
    variables$states$.update()

    # event clocks tick
    for (event in events[events_2_loop]) {
      event$.tick()
    }
    for (event in events$scheduled_dose) {
      event$.tick()
    }

    # print progress?
    if (progress) {
      setTxtProgressBar(pb = pb,value = t)
    }
  }
  # end loop, close out pb
  if (progress) {
    close(con = pb)
  }
}
