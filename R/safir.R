#' safir: squire and friends individual rewrite
#'
#' @description
#' safir is an individual-based simulation model of COVD-19 based on squire.
#'
#' @docType package
#' @name safir
#'
#' @importFrom stats setNames
#' @importFrom utils getFromNamespace
#'
"_PACKAGE"

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @useDynLib safir
## usethis namespace: end
NULL

execute_any_process <- utils::getFromNamespace("execute_any_process", "individual")
execute_process <- utils::getFromNamespace("execute_process", "individual")
nimue_parameters <- utils::getFromNamespace("parameters", "nimue")
nimue_probs <- utils::getFromNamespace("probs", "nimue")
nimue_durs <- utils::getFromNamespace("durs", "nimue")
nimue_vaccine_pars <- utils::getFromNamespace("vaccine_pars", "nimue")
nimue_eligable_for_second <- utils::getFromNamespace("eligable_for_second", "nimue")
nimue_target_pop <- utils::getFromNamespace("target_pop", "nimue")
