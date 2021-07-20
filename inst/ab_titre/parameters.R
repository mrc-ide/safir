mu_ab_d1 <- 0.14 # mean titre dose 1
mu_ab_d2 <- 2.37 # mean titre dose 2
ab_50 <- 0.2 # titre relative to convalescent required to provide 50% protection from infection, on linear scale
ab_50_severe <- 0.03
std10 <- 0.44 # Pooled standard deviation of antibody level on log10 data
k <- 2.94 # shape parameter of efficacy curve
max_t <- 730 # number of days to model
t_d2 <- 84 # timing of second dose relative to first
hl_s <- 108 # Half life of antibody decay - short
hl_l <- 3650 # Half life of antibody decay - long
period_s <- 250
t_period_l <- 365 # Time point at which to have switched to longest half-life
dr_s <- -log(2)/hl_s # Corresponding decay rate in days for half life above
dr_l <- -log(2)/hl_l

nt <- NULL
t <- 0:max_t
time_to_decay <- t_period_l - period_s # time in days to reach longest half-life

# vector of decay rates over time: first we have a period of fast decay, then gradually shift to period of long decay
dr_vec <- c(rep(dr_s, period_s),
            seq(dr_s, dr_l, length.out = time_to_decay),
            rep(dr_l, (length(t) - t_period_l)))



hl_s <- 108 # Half life of antibody decay - short
hl_l <- 3650 # Half life of antibody decay - long
period_s <- 250
t_period_l <- 365 # Time point at which to have switched to longest half-life
time_to_decay <- t_period_l - period_s # time in days to reach longest half-life
dr_s <- -log(2)/hl_s # Corresponding decay rate in days for half life above
dr_l <- -log(2)/hl_l

dr_vec <- c(rep(dr_s, period_s),
            seq(dr_s, dr_l, length.out = time_to_decay),
            dr_l)

parameters <- list(
  dr_vec = dr_vec,
  mu_ab = c(mu_ab_d1, mu_ab_d2)
)


