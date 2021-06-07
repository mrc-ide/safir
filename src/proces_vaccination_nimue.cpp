/* --------------------------------------------------------------------------------
 *  vaccination process for nimue vaccine model
 *  Sean L. Wu (slwood89@gmail.com)
 *  June 2021
 -------------------------------------------------------------------------------- */

#include <Rcpp.h>
#include <individual.h>
#include "../inst/include/utils.hpp"

//' @title C++ infection process (nimue vaccine model)
//' @description this is an internal function, you should use the R interface
//' for type checking, \code{\link{infection_process_cpp}}
//' @param parameters a list of parameters from \code{\link{get_parameters_nimue}}
//' @param states a \code{\link[individual]{CategoricalVariable}}
//' @param eligible a \code{\link[individual]{Bitset}}
//' @param vaccinated a \code{\link[individual]{Bitset}}
//' @param empty a \code{\link[individual]{Bitset}}
//' @param discrete_age a \code{\link[individual]{IntegerVariable}}
//' @param v0_to_v1v2 a \code{\link[individual]{TargetedEvent}}
//' @param dt size of time step
// [[Rcpp::export]]
Rcpp::XPtr<process_t> vaccination_process_nimue_cpp_internal(
    Rcpp::List parameters,
    Rcpp::XPtr<CategoricalVariable> states,
    Rcpp::XPtr<individual_index_t> eligible,
    Rcpp::XPtr<individual_index_t> vaccinated,
    Rcpp::XPtr<individual_index_t> empty,
    Rcpp::XPtr<IntegerVariable> discrete_age,
    Rcpp::XPtr<TargetedEvent> v0_to_v1v2,
    const double dt
) {

  // states that may receive vaccination
  std::vector<std::string> vaxx_states = {"S", "E", "R"};

  // fixed parameters
  int N_age = Rcpp::as<int>(parameters["N_age"]);
  int N_prioritisation_steps = Rcpp::as<int>(parameters["N_prioritisation_steps"]);
  Rcpp::NumericMatrix vaccine_coverage_mat = Rcpp::as<Rcpp::NumericMatrix>(parameters["vaccine_coverage_mat"]);

  // make these once, modify them by specifying mutable value captures (can't use
  // reference captures since their lifetime ends when the func returns)
  std::vector<double> pr(N_age, 0.);
  Rcpp::IntegerMatrix vaccination_target_mat(N_prioritisation_steps , N_age);
  std::vector<int> vaccine_target_vec(N_prioritisation_steps, 0);
  std::vector<int> vaccination_target(N_age, 0);

  // the process lambda
  return Rcpp::XPtr<process_t>(
    new process_t([parameters, states, eligible, vaccinated, empty, discrete_age, v0_to_v1v2, dt, vaxx_states, N_age, N_prioritisation_steps, vaccine_coverage_mat, pr, vaccination_target_mat, vaccine_target_vec, vaccination_target](size_t t) mutable {

      // current day (subtract one for zero-based indexing)
      size_t tnow = std::ceil((double)t * dt) - 1.;

      double mv = get_vector_cpp(parameters["vaccine_set"], tnow);

      if (mv > 0.) {

        // calculate prioritisation step and which age groups are eligible right now
        for (auto i = 0u; i < N_age; ++i) {
          pr[i] = get_proportion_vaccinated_nimue_internal(discrete_age, vaccinated, i+1);
        }

        for (auto p = 0u; p < N_prioritisation_steps; ++p) {
          for (auto i = 0; i < N_age; ++i) {
            vaccination_target_mat(p, i) = (pr[i] < vaccine_coverage_mat(p, i))? 1 : 0;
          }
        }

        // an entire row summing to zero means that step has been completed
        for (auto p = 0u; p < N_prioritisation_steps; ++p) {
          int rowsum{0};
          for (auto i = 0; i < N_age; ++i) {
            rowsum += vaccination_target_mat(p, i);
          }
          vaccine_target_vec[p] = (rowsum == 0)? 1 : 0;
        }
        int current_index = std::min(std::accumulate(vaccine_target_vec.begin(), vaccine_target_vec.end(), 0) + 1, N_prioritisation_steps);
        current_index -= 1; // 0-based indexing

        for (auto i = 0u; i < N_age; ++i) {
          vaccination_target[i] = vaccination_target_mat(current_index, i);
        }

        // if no vaccination targets remain don't run the code to distribute vaccines
        if (std::accumulate(vaccination_target.begin(), vaccination_target.end(), 0) > 0) {

          Rcpp::Rcout << "distributing vaxx on step: " << current_index + 1 << " --- \n";

          // clear out eligible
          *eligible &= *empty;

          // get SER in eligible ages
          individual_index_t SER = states->get_index_of(vaxx_states);
          std::vector<int> target_ages;
          for (auto i = 0u; i < N_age; ++i) {
            if (vaccination_target[i] > 0) {
              target_ages.emplace_back(i+1);
            }
          }
          individual_index_t target_ages_bitset = discrete_age->get_index_of_set(target_ages);
          SER &= target_ages_bitset;

          // set who is eligible: SER people in an age group in this priority step AND unvaccinated
          *eligible |= SER;
          *eligible &= ~(*vaccinated);

          // calc rate of vaccination now
          double vr_den = eligible->size();
          double vr = std::min(mv/vr_den, 1.);

          // sample who gets vaccinated
          bitset_sample_internal(*eligible, Rf_pexp(vr * dt, 1., 1, 0));
          if (eligible->size() > 0) {
            v0_to_v1v2->schedule(*eligible, 0.);
          }

        } // end check

      } // end if

    }),
    true
  );
};
