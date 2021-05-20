# --------------------------------------------------------------------------------
#   
# --------------------------------------------------------------------------------

#' @title Schedule events for individuals at initialisation
#' @export
setup_events <- function(
    parameters,
    events,
    variables,
    dt,
    vaccines = NULL
) {

    # E individuals
    bset_E <- variables$states$get_index_of("E")
    if (bset_E$size() > 0) {
        init_fn <- create_exposure_scheduler(
            events = events,
            variables = variables,
            parameters = parameters,
            dt = dt,
            vaccines = vaccines
        )
        init_fn(timestep = 1, to_move = bset_E)
    }

    # IMild
    bset_IMild <- variables$states$get_index_of("IMild")
    if (bset_IMild$size() > 0) {
        init_fn <- create_event_scheduler_listener(
            event = events$recovery,
            duration = parameters$dur_IMild,
            func = make_rexp_simple,
            shift = 1L,
            dt = dt
        )
        init_fn(timestep = 1, to_move = bset_IMild)
    }

    # IAsymp
    bset_IAsymp <- variables$states$get_index_of("IAsymp")
    if (bset_IAsymp$size() > 0) {
        init_fn <- create_event_scheduler_listener(
            event = events$recovery,
            duration = parameters$dur_IAsymp,
            func = make_rexp_simple,
            shift = 1L,
            dt = dt
        )
        init_fn(timestep = 1, to_move = bset_IAsymp)
    }

    # ICase
    bset_ICase <- variables$states$get_index_of("ICase")
    if (bset_ICase$size() > 0) {
        init_fn <- create_event_scheduler_listener(
            event = events$hospitilisation,
            duration = parameters$dur_ICase,
            func = make_rerlang,
            shift = 1L,
            dt = dt
        )
        init_fn(timestep = 1, to_move = bset_ICase)
    }

    # IMVGetLive
    bset_IMVGetLive <- variables$states$get_index_of("IMVGetLive")
    if (bset_IMVGetLive$size() > 0) {
        init_fn <- create_event_scheduler_listener(
            event = events$stepdown,
            duration = parameters$dur_get_mv_survive,
            func = make_rerlang,
            shift = 1L,
            dt = dt
        )
        init_fn(timestep = 1, to_move = bset_IMVGetLive)
    }

    # IMVGetDie
    bset_IMVGetDie <- variables$states$get_index_of("IMVGetDie")
    if (bset_IMVGetDie$size() > 0) {
        init_fn <- create_event_scheduler_listener(
            event = events$death,
            duration = parameters$dur_get_mv_die,
            func = make_rerlang,
            shift = 1L,
            dt = dt
        )
        init_fn(timestep = 1, to_move = bset_IMVGetDie)
    }

    # IMVNotGetLive
    bset_IMVNotGetLive <- variables$states$get_index_of("IMVNotGetLive")
    if (bset_IMVNotGetLive$size() > 0) {
        init_fn <- create_event_scheduler_listener(
            event = events$recovery,
            duration = parameters$dur_not_get_mv_survive,
            func = make_rerlang,
            shift = 1L,
            dt = dt
        )
        init_fn(timestep = 1, to_move = bset_IMVNotGetLive)
    }

    # IMVNotGetDie
    bset_IMVNotGetDie <- variables$states$get_index_of("IMVNotGetDie")
    if (bset_IMVNotGetDie$size() > 0) {
        init_fn <- create_event_scheduler_listener(
            event = events$death,
            duration = parameters$dur_not_get_mv_die,
            func = make_rerlang,
            shift = 1L,
            dt = dt
        )
        init_fn(timestep = 1, to_move = bset_IMVNotGetDie)
    }

    # IOxGetLive
    bset_IOxGetLive <- variables$states$get_index_of("IOxGetLive")
    if (bset_IOxGetLive$size() > 0) {
        init_fn <- create_event_scheduler_listener(
            event = events$recovery,
            duration = parameters$dur_get_ox_survive,
            func = make_rerlang,
            shift = 1L,
            dt = dt
        )
        init_fn(timestep = 1, to_move = bset_IOxGetLive)
    }

    # IOxGetDie
    bset_IOxGetDie <- variables$states$get_index_of("IOxGetDie")
    if (bset_IOxGetDie$size() > 0) {
        init_fn <- create_event_scheduler_listener(
            event = events$death,
            duration = parameters$dur_get_ox_die,
            func = make_rerlang,
            shift = 1L,
            dt = dt
        )
        init_fn(timestep = 1, to_move = bset_IOxGetDie)
    }

    # IOxNotGetLive
    bset_IOxNotGetLive <- variables$states$get_index_of("IOxNotGetLive")
    if (bset_IOxNotGetLive$size() > 0) {
        init_fn <- create_event_scheduler_listener(
            event = events$recovery,
            duration = parameters$dur_not_get_ox_survive,
            func = make_rerlang,
            shift = 1L,
            dt = dt
        )
        init_fn(timestep = 1, to_move = bset_IOxNotGetLive)
    }

    # IOxNotGetDie
    bset_IOxNotGetDie <- variables$states$get_index_of("IOxNotGetDie")
    if (bset_IOxNotGetDie$size() > 0) {
        init_fn <- create_event_scheduler_listener(
            event = events$death,
            duration = parameters$dur_not_get_ox_die,
            func = make_rerlang,
            shift = 1L,
            dt = dt
        )
        init_fn(timestep = 1, to_move = bset_IOxNotGetDie)
    }

    # R
    if (is.finite(parameters$dur_R)) {
        bset_R <- variables$states$get_index_of("R")
        if (bset_R$size() > 0){
            init_fn <- create_event_scheduler_listener(
                event = events$immunity_loss,
                duration = parameters$dur_R,
                func = make_rerlang,
                shift = 1L,
                dt = dt
            )
            init_fn(timestep = 1, to_move = bset_R)
        }
    }

    # IRec
    bset_IRec <- variables$states$get_index_of("IRec")
    if (bset_IRec$size() > 0) {
        init_fn <- create_event_scheduler_listener(
            event = events$recovery,
            duration = parameters$dur_rec,
            func = make_rerlang,
            shift = 1L,
            dt = dt
        )
        init_fn(timestep = 1, to_move = bset_IRec)
    }

}