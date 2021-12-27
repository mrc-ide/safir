# renderer objects for hospitalization and ICU outcomes

#' @title Create list of renderer objects to track hospitalization and ICU incidence
#' @param parameters named list of model parameters
#' @return a named list of [individual::Render] objects
#' @importFrom individual Render
#' @export
create_hosp_renderers <- function(parameters) {

  time_steps <- parameters$time_steps
  stopifnot(!is.null(time_steps))

  list(
    "ICU_get_live" = Render$new(timesteps = time_steps),
    "ICU_get_die" = Render$new(timesteps = time_steps),
    "hosp_get_live" = Render$new(timesteps = time_steps),
    "hosp_get_die" = Render$new(timesteps = time_steps),
    "ICU_not_get_live" = Render$new(timesteps = time_steps),
    "ICU_not_get_die" = Render$new(timesteps = time_steps),
    "hosp_not_get_live" = Render$new(timesteps = time_steps),
    "hosp_not_get_die" = Render$new(timesteps = time_steps)
  )
}


#' @title Attach listeners to track hospitalization and ICU incidence
#' @param renderers a named list of [individual::Render] objects
#' @param events a named list of [individual::TargetedEvent] objects
#' @export
attach_hosp_listeners <- function(renderers, events) {

  event_tags <- c("imv_get_live", "imv_get_die", "iox_get_live", "iox_get_die", "imv_not_get_live", "imv_not_get_die", "iox_not_get_live", "iox_not_get_die")
  stopifnot(vapply(X = events[event_tags], FUN = function(x) {inherits(x, "TargetedEvent")}, FUN.VALUE = logical(1)))
  stopifnot(vapply(X = renderers, FUN = function(x) {inherits(x, "Render")}, FUN.VALUE = logical(1)))

  # need ICU, gets bed, lives
  events$imv_get_live$add_listener(
    function(timestep, target) {
      renderers$ICU_get_live$render(name = "ICU_get_live", value = target$size(), timestep = timestep)
    }
  )

  # need ICU, gets bed, dies
  events$imv_get_die$add_listener(
    function(timestep, target) {
      renderers$ICU_get_die$render(name = "ICU_get_die", value = target$size(), timestep = timestep)
    }
  )

  # need hosp, gets bed, lives
  events$iox_get_live$add_listener(
    function(timestep, target) {
      renderers$hosp_get_live$render(name = "hosp_get_live", value = target$size(), timestep = timestep)
    }
  )

  # need hosp, gets bed, dies
  events$iox_get_die$add_listener(
    function(timestep, target) {
      renderers$hosp_get_die$render(name = "hosp_get_die", value = target$size(), timestep = timestep)
    }
  )

  # need ICU, doesn't get bed, lives
  events$imv_not_get_live$add_listener(
    function(timestep, target) {
      renderers$ICU_not_get_live$render(name = "ICU_not_get_live", value = target$size(), timestep = timestep)
    }
  )

  # need ICU, doesn't get bed, dies
  events$imv_not_get_die$add_listener(
    function(timestep, target) {
      renderers$ICU_not_get_die$render(name = "ICU_not_get_die", value = target$size(), timestep = timestep)
    }
  )

  # need hosp, doesn't get bed, lives
  events$iox_not_get_live$add_listener(
    function(timestep, target) {
      renderers$hosp_not_get_live$render(name = "hosp_not_get_live", value = target$size(), timestep = timestep)
    }
  )

  # need hosp, doesn't get bed, dies
  events$iox_not_get_die$add_listener(
    function(timestep, target) {
      renderers$hosp_not_get_die$render(name = "hosp_not_get_die", value = target$size(), timestep = timestep)
    }
  )
}


#' @title Process rendering output of hospitalization and ICU incidence
#' @description Process output from each rendering object and return a single
#' [data.table::data.table] object, where incidence has been aggregated by
#' day.
#' @param renderers a named list of [individual::Render] objects from [safir::create_hosp_renderers]
#' @param parameters a named list of model parameters
#' @importFrom data.table as.data.table .SD setnames
#' @export
process_hosp_renderers <- function(renderers, parameters) {
  # process hospital renderers
  out <- mapply(FUN = function(x, name) {
    df <- x$to_dataframe()
    if (ncol(df) == 1L) {
      # attach zero col with right name
      df[[name]] = rep(0, nrow(df))
    }
    df[which(is.na(df[, 2])), 2] <- 0
    return(df)
  }, x = renderers, name = names(renderers), SIMPLIFY = FALSE, USE.NAMES = FALSE)
  out <- do.call(cbind, out)
  out <- out[, -which(names(out) == "timestep")[-1]]
  # go to days
  out$timestep <- as.integer(floor(out$timestep * parameters$dt) + 1)
  out <- as.data.table(out)
  out <- out[, lapply(.SD, sum) , by = "timestep", .SDcols = !c("timestep")]
  setnames(x = out, old = "timestep", new = "day")
  return(out)
}
