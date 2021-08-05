# --------------------------------------------------------------------------------
#   infection process for vaccination model (multiple doses, no types)
#   Sean L. Wu (slwood89@gmail.com)
#   July 2021
# --------------------------------------------------------------------------------


#' @title Get vaccine efficacy and Ab titre parameters
#' @param vaccine which vaccine? should be one of: "Pfizer", "AstraZeneca", "Sinovac", "Moderna"
#' @param max_dose maximum number of doses
#' @param correlated are doses correlated?
#' @description Get parameters for vaccine efficacy and antibody titre decay rate.
#' @export
get_vaccine_ab_titre_parameters <- function(vaccine, max_dose = 2, correlated = FALSE) {
  stopifnot(max_dose == 2)
  stopifnot(is.logical(correlated))
  stopifnot(vaccine %in% c("Pfizer", "AstraZeneca", "Sinovac", "Moderna"))

  mu_ab_list <- data.frame(name = c("Pfizer", "AstraZeneca", "Sinovac", "Moderna"),
                           mu_ab_d1 = c(13/94, 1/59, 28/164, ((185+273)/2)/321),
                           mu_ab_d2 = c(223/94, 32/59, 28/164, 654/158))

  hl_s <- 108 # Half life of antibody decay - short
  hl_l <- 3650 # Half life of antibody decay - long
  period_s <- 250
  t_period_l <- 365 # Time point at which to have switched to longest half-life
  time_to_decay <- t_period_l - period_s # time in days to reach longest half-life
  dr_s <- -log(2)/hl_s # Corresponding decay rate in days for half life above
  dr_l <- -log(2)/hl_l

  mu_ab_d1 <- mu_ab_list[mu_ab_list$name == vaccine, "mu_ab_d1"] # mean titre dose 1
  mu_ab_d2 <- mu_ab_list[mu_ab_list$name == vaccine, "mu_ab_d2"] # mean titre dose 2
  ab_50 <- 0.2 # titre relative to convalescent required to provide 50% protection from infection, on linear scale
  ab_50_severe <- 0.03
  std10 <- 0.44 # Pooled standard deviation of antibody level on log10 data
  k <- 2.94 # shape parameter of efficacy curve

  dr_vec <- c(
    rep(dr_s, period_s),
    seq(dr_s, dr_l, length.out = time_to_decay),
    dr_l
  )

  parameters <- list(
    dr_vec = dr_vec,
    mu_ab = c(mu_ab_d1, mu_ab_d2),
    std10 = std10,
    ab_50 = ab_50,
    ab_50_severe = ab_50_severe,
    k = k,
    correlated = correlated
  )
  return(parameters)
}



