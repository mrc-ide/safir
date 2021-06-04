# --------------------------------------------------
#   Categorical infection states (squire)
#   Sean L. Wu (slwood89@gmail.com)
#   May 2021
# --------------------------------------------------


#' @title Define model states
#' @description
#' Create_states creates and initialises the human states for the model
#'
#' @title Create and initialise states
#'
#' @param psq the model parameters
#'
#' @return list of states
#' @export
create_state_variables <- function(psq) {

  # Sum up the states from squire
  initial_counts <- c(
    S = sum(psq$S_0),
    E = sum(psq$E1_0, psq$E2_0),
    IMild = sum(psq$IMild_0),
    IAsymp = sum(psq$IAsymp_0),
    ICase = sum(psq$ICase1_0, psq$ICase2_0),
    IOxGetLive = sum(psq$IOxGetLive1_0, psq$IOxGetLive2_0),
    IOxGetDie = sum(psq$IOxGetDie1_0, psq$IOxGetDie2_0),
    IOxNotGetLive = sum(psq$IOxNotGetLive1_0, psq$IOxNotGetLive2_0),
    IOxNotGetDie = sum(psq$IOxNotGetDie1_0, psq$IOxNotGetDie2_0),
    IMVGetLive = sum(psq$IMVGetLive1_0, psq$IMVGetLive2_0),
    IMVGetDie = sum(psq$IMVGetDie1_0, psq$IMVGetDie2_0),
    IMVNotGetLive = sum(psq$IMVNotGetLive1_0, psq$IMVNotGetLive2_0),
    IMVNotGetDie = sum(psq$IMVNotGetDie1_0, psq$IMVNotGetDie2_0),
    IRec = sum(psq$IRec1_0, psq$IRec2_0),
    R = sum(psq$R_0),
    D = sum(psq$D_0)
  )

  # Define state variables
  states <- individual::CategoricalVariable$new(
    names(initial_counts),
    rep(names(initial_counts), times = as.integer(initial_counts))
  )

  return(states)
}
