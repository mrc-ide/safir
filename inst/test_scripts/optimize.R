


schedule_outcome <- function(
  target,
  prob_successful,
  success_event,
  failure_event
) {
  success <- bernoulli_multi_p(prob_successful)

  if(sum(success) > 0) {
    to_success <- individual::filter_bitset(target, which(success))
    success_event$schedule(to_success, delay = 0)
  }

  if(sum(!success) > 0) {
    to_failure <- individual::filter_bitset(target, which(!success))
    failure_event$schedule(to_failure, delay = 0)
  }
}

library(microbenchmark)
library(ggplot2)

# prob_successful <- runif(n = 1e4)
prob_successful <- c(rep(1,5e3),rep(0,5e3))
target_0 <- Bitset$new(size = 1e5)$insert(sample.int(n = 1e5,size = 1e4,replace = FALSE))

# original
target <- target_0$copy()
success <- safir:::bernoulli_multi_p(prob_successful)

if(sum(success) > 0) {
  to_success <- individual::filter_bitset(target, which(success))
}

if(sum(!success) > 0) {
  to_failure <- individual::filter_bitset(target, which(!success))
}

# new
target <- target_0$copy()
success <- target$copy()

success$sample(rate = prob_successful)
failure <- target$copy()$set_difference(success)


xx <- microbenchmark(
  "original" = {target <- target_0$copy()
  success <- safir:::bernoulli_multi_p(prob_successful)

  if(sum(success) > 0) {
    to_success <- individual::filter_bitset(target, which(success))
  }

  if(sum(!success) > 0) {
    to_failure <- individual::filter_bitset(target, which(!success))
  }},
  "new" = {
    target <- target_0$copy()
  success <- target$copy()

  success$sample(rate = prob_successful)
  failure <- target$copy()$set_difference(success)}
)
