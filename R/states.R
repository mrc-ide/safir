#' @title Define model states
#' @description
#' Create_states creates and initialises the human states for the model
#'
#' @title Create and initialise states
#'
#' @param psq the model parameters
#'
#' @noRd
#' @return list of states
create_states <- function(psq) {

  # Sum up the states from squire
  S <- sum(psq$S_0)
  E <- sum(psq$E1_0, psq$E2_0)
  IMild <- sum(psq$IMild_0)
  IAsymp <- sum(psq$IAsymp_0)
  ICase <- sum(psq$ICase1_0, psq$ICase2_0)
  IOxGetLive <- sum(psq$IOxGetLive1_0, psq$IOxGetLive2_0)
  IOxGetDie <- sum(psq$IOxGetDie1_0, psq$IOxGetDie2_0)
  IOxNotGetLive <- sum(psq$IOxNotGetLive1_0, psq$IOxNotGetLive2_0)
  IOxNotGetDie <- sum(psq$IOxNotGetDie1_0, psq$IOxNotGetDie2_0)
  IMVGetLive <- sum(psq$IMVGetLive1_0, psq$IMVGetLive2_0)
  IMVGetDie <- sum(psq$IMVGetDie1_0, psq$IMVGetDie2_0)
  IMVNotGetLive <- sum(psq$IMVNotGetLive1_0, psq$IMVNotGetLive2_0)
  IMVNotGetDie <- sum(psq$IMVNotGetDie1_0, psq$IMVNotGetDie2_0)
  IRec <- sum(psq$IRec1_0, psq$IRec2_0)
  R <- sum(psq$R_0)
  D <- sum(psq$D_0)

  # Create states for Humans
  list(
    S = individual::State$new("S", S),
    E = individual::State$new("E", E),
    IMild = individual::State$new("IMild", IMild),
    IAsymp = individual::State$new("IAsymp", IAsymp),
    ICase = individual::State$new("ICase", ICase),
    IOxGetLive = individual::State$new("IOxGetLive", IOxGetLive),
    IOxGetDie = individual::State$new("IOxGetDie", IOxGetDie),
    IOxNotGetLive = individual::State$new("IOxNotGetLive", IOxNotGetLive),
    IOxNotGetDie = individual::State$new("IOxNotGetDie", IOxNotGetDie),
    IMVGetLive = individual::State$new("IMVGetLive", IMVGetLive),
    IMVGetDie = individual::State$new("IMVGetDie", IMVGetDie),
    IMVNotGetLive = individual::State$new("IMVNotGetLive", IMVNotGetLive),
    IMVNotGetDie = individual::State$new("IMVNotGetDie", IMVNotGetDie),
    IRec = individual::State$new("IRec", IRec),
    R = individual::State$new("R", R),
    D = individual::State$new("D", D)
  )
}
