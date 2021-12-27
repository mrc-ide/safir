#include "../inst/include/events_hospital.hpp"

// return a lambda that gets proabilities given (time, ages)
get_probs_fn make_get_probs(SEXP probs) {
  // for time varying (daily) probs
  if (Rf_isMatrix(probs)) {
    int n_ages = Rf_nrows(probs);
    double* probs_ptr = REAL(probs);

    return [n_ages, probs_ptr] (const size_t day, const std::vector<int>& ages) -> std::vector<double> {
      std::vector<double> probs(ages.size());
      for (auto i = 0u; i < ages.size(); ++i) {
        probs[i] = probs_ptr[day * n_ages + ages[i]];
      }
      return probs;
    };
  } else {
    // probs do not vary
    double* probs_ptr = REAL(probs);
    return [probs_ptr] (const size_t day, const std::vector<int>& ages) -> std::vector<double> {
      std::vector<double> probs(ages.size());
      for (auto i = 0u; i < ages.size(); ++i) {
        probs[i] = probs_ptr[ages[i]];
      }
      return probs;
    };
  }
};


// make a lambda that allocates persons to treatment and returns a Bitset
allocate_treatment_fn get_allocate_treatment(
    Rcpp::XPtr<CategoricalVariable> states,
    const std::vector<std::string> treated_states
) {
  return [states, treated_states] (const individual_index_t& need_treatment, const int limit) -> individual_index_t {
    // allocation
    individual_index_t occupied = states->get_index_of(treated_states);
    int available = limit - static_cast<int>(occupied.size());

    // who is getting treatment from available
    if (need_treatment.size() <= available) {
      return need_treatment;
    }

    // allocate with scarcity
    int k = std::max(0, available);
    if (k > 0) {
      // k persons get treatment
      individual_index_t get_treatment(need_treatment);
      bitset_choose_internal(get_treatment, k);
      return get_treatment;
    } else {
      // nobody gets treatment
      return individual_index_t(need_treatment.max_size());
    }
  };
};


// make a lambda that schedules individuals for success or failure (Bernoulli event)
schedule_outcome_fn get_schedule_outcome(
    Rcpp::XPtr<TargetedEvent> success_event,
    Rcpp::XPtr<TargetedEvent> failure_event
) {
  return [success_event, failure_event] (const individual_index_t& target, const std::vector<double>& prob_successful) -> void {
    individual_index_t success(target);
    bitset_sample_multi_internal(success, prob_successful.begin(), prob_successful.end());

    individual_index_t failure(target);
    failure &= !success; // set difference with success

    if (success.size() > 0) {
      success_event->schedule(success, 0.0);
    }

    if (failure.size() > 0) {
      failure_event->schedule(failure, 0.0);
    }

  };
};


// // [[Rcpp::export]]
// Rcpp::XPtr<targeted_listener_t> create_hospital_scheduler_listener_cpp_internal(
//     Rcpp::List parameters,
//     Rcpp::XPtr<CategoricalVariable> states,
//     Rcpp::XPtr<IntegerVariable> discrete_age,
//     Rcpp::XPtr<TargetedEvent> imv_get_die,
//     Rcpp::XPtr<TargetedEvent> imv_get_live,
//     Rcpp::XPtr<TargetedEvent> imv_not_get_die,
//     Rcpp::XPtr<TargetedEvent> imv_not_get_live,
//     Rcpp::XPtr<TargetedEvent> iox_get_die,
//     Rcpp::XPtr<TargetedEvent> iox_get_live,
//     Rcpp::XPtr<TargetedEvent> iox_not_get_die,
//     Rcpp::XPtr<TargetedEvent> iox_not_get_live
// ) {
//   double dt = Rcpp::as<double>(parameters["dt"]);
//   int ICU_beds = Rcpp::as<int>(parameters["ICU_beds"]);
//   int hosp_beds = Rcpp::as<int>(parameters["hosp_beds"]);
//
//   // ICU states
//   std::vector<std::string> ICU_states = {"IMVGetDie", "IMVGetLive"};
//
//   // hosp states
//   std::vector<std::string> hosp_states = {"IOxGetDie", "IOxGetLive", "IRec"};
//
//   // get probabilities
//   get_probs_fn get_prob_severe = make_get_probs(parameters["prob_severe"]);
//   get_probs_fn get_prob_severe_death_treatment = make_get_probs(parameters["prob_severe_death_treatment"]);
//   get_probs_fn get_prob_severe_death_no_treatment = make_get_probs(parameters["prob_severe_death_no_treatment"]);
//   get_probs_fn get_prob_non_severe_death_treatment = make_get_probs(parameters["prob_non_severe_death_treatment"]);
//   get_probs_fn get_prob_non_severe_death_no_treatment = make_get_probs(parameters["prob_non_severe_death_no_treatment"]);
//
//   // allocate treatment
//   allocate_treatment_fn allocate_ICU = get_allocate_treatment(states, ICU_states);
//   allocate_treatment_fn allocate_hosp = get_allocate_treatment(states, hosp_states);
//
//   // schedule outcomes
//   schedule_outcome_fn schedule_outcome_ICU_get = get_schedule_outcome(imv_get_die, imv_get_live);
//   schedule_outcome_fn schedule_outcome_ICU_not_get = get_schedule_outcome(imv_not_get_die, imv_not_get_live);
//   schedule_outcome_fn schedule_outcome_hosp_get = get_schedule_outcome(iox_get_die, iox_get_live);
//   schedule_outcome_fn schedule_outcome_hosp_not_get = get_schedule_outcome(iox_not_get_die, iox_not_get_live);
//
//   // generate the lambda fn
//   return Rcpp::XPtr<targeted_listener_t>(
//     new targeted_listener_t([
//                               dt, ICU_beds, hosp_beds,
//                               states, discrete_age,
//                               get_prob_severe, get_prob_severe_death_treatment, get_prob_severe_death_no_treatment, get_prob_non_severe_death_treatment, get_prob_non_severe_death_no_treatment,
//                               allocate_ICU, allocate_hosp,
//                               schedule_outcome_ICU_get, schedule_outcome_ICU_not_get, schedule_outcome_hosp_get, schedule_outcome_hosp_not_get
//                               ](size_t timestep, const individual_index_t& hospitalised) -> void {
//
//       // current day (subtract one for zero-based indexing)
//       size_t day = std::ceil((double)timestep * dt) - 1.;
//
//       std::vector<int> disc_ages = discrete_age->get_values(hospitalised);
//       std::vector<double> prob_severe = get_prob_severe(day, disc_ages);
//
//       individual_index_t need_mv(hospitalised);
//       bitset_sample_multi_internal(need_mv, prob_severe.begin(), prob_severe.end());
//
//       individual_index_t need_ox(hospitalised);
//       need_ox &= !need_mv;
//
//       // individuals requiring ICU
//       if (need_mv.size() > 0) {
//
//         // number of people who can get ICU bed
//         individual_index_t mv_get = allocate_ICU(need_mv, ICU_beds);
//
//         // schedule for those getting a ICU bed
//         if (mv_get.size() > 0) {
//           std::vector<int> mv_get_ages = discrete_age->get_values(mv_get);
//           std::vector<double> prob_death = get_prob_severe_death_treatment(day, mv_get_ages);
//           schedule_outcome_ICU_get(mv_get, prob_death);
//         }
//
//         // schedule for those not getting a ICU bed
//         need_mv &= !mv_get;
//         if (need_mv.size() > 0) {
//           std::vector<int> mv_not_get_ages = discrete_age->get_values(need_mv);
//           std::vector<double> prob_death = get_prob_severe_death_no_treatment(day, mv_not_get_ages);
//           schedule_outcome_ICU_not_get(need_mv, prob_death);
//         }
//
//       }
//
//       // individuals requiring hosp
//       if (need_ox.size() > 0) {
//
//         // number of people who can get hosp bed
//         individual_index_t ox_get = allocate_hosp(need_ox, hosp_beds);
//
//         // schedule for those getting hosp bed
//         if (ox_get.size() > 0) {
//           std::vector<int> ox_get_ages = discrete_age->get_values(ox_get);
//           std::vector<double> prob_death = get_prob_non_severe_death_treatment(day, ox_get_ages);
//           schedule_outcome_hosp_get(ox_get, prob_death);
//         }
//
//         // schedule for those not getting hosp bed
//         need_ox &= !ox_get;
//         if (need_ox.size() > 0) {
//           std::vector<int> ox_not_get_ages = discrete_age->get_values(need_ox);
//           std::vector<double> prob_death = get_prob_non_severe_death_no_treatment(day, ox_not_get_ages);
//           schedule_outcome_hosp_not_get(need_ox, prob_death);
//         }
//       }
//
//     }), // end of lambda
//     true
//   );
// };
