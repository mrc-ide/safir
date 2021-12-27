// schedule events upon needing hospitalization

#ifndef EVENTS_HOSP_HPP
#define EVENTS_HOSP_HPP

#include <vector>
#include <functional>
#include <string>
#include <Rcpp.h>
#include <individual.h>

// function that returns a vector of probabilities (args: timestep, ages)
using get_probs_fn = std::function<std::vector<double>(const size_t, const std::vector<int>&)>;

get_probs_fn make_get_probs(SEXP probs);

// function that allocates treatment (args: need_treatment, limit)
using allocate_treatment_fn = std::function<individual_index_t(const individual_index_t&, const int)>;

allocate_treatment_fn get_allocate_treatment(
    Rcpp::XPtr<CategoricalVariable> states,
    const std::vector<std::string> treated_states
);

// function to schedule outcomes (args: target, prob_successful)
using schedule_outcome_fn = std::function<void(const individual_index_t&, const std::vector<double>&)>;

schedule_outcome_fn get_schedule_outcome(
  Rcpp::XPtr<TargetedEvent> success_event,
  Rcpp::XPtr<TargetedEvent> failure_event
);

#endif
