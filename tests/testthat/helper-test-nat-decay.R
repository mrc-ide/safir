draw_nt <- function(parameters, n, tmax) {

  # get pars out
  hl_s <- parameters$hl_s # Half life of antibody decay - short
  hl_l <- parameters$hl_l # Half life of antibody decay - long
  period_s <- parameters$period_s
  t_period_l <- parameters$t_period_l # Time point at which to have switched to longest half-life

  time_to_decay <- t_period_l - period_s # time in days to reach longest half-life
  dr_s <- -log(2)/hl_s # Corresponding decay rate in days for half life above
  dr_l <- -log(2)/hl_l

  ab_50 <- parameters$ab_50 # titre relative to convalescent required to provide 50% protection from infection, on linear scale
  ab_50_severe <- parameters$ab_50_severe
  k <- parameters$k # shape parameter of efficacy curve

  mu_ab_d1 <- parameters$mu_ab[1] # mean titre dose 1
  std10 <- parameters$std10 # Pooled standard deviation of antibody level on log10 data

  t <- 0:tmax

  # vector of decay rates over time: first we have a period of fast decay, then gradually shift to period of long decay
  dr_vec <- c(rep(dr_s, period_s),
              seq(dr_s, dr_l, length.out = time_to_decay),
              rep(dr_l, (length(t) - t_period_l)))

  z1 <- rnorm(n, log10(mu_ab_d1), std10)

  # initiate titre vector
  nt <- matrix(data = NaN,nrow = length(t),ncol = n)
  nt[1, ] <- log(10^z1)

  # decay antibodies over time on natural log scale
  for (i in (2:length(t))){
    nt[i, ] <- nt[i-1, ] + dr_vec[i-1]
  }

  return(list(nt = nt, z1 = z1))
}
