/* --------------------------------------------------------------------------------
 *  infection process for squire transmission model
 *  Sean L. Wu (slwood89@gmail.com)
 *  June 2021
 -------------------------------------------------------------------------------- */

#include <Rcpp.h>
#include <individual.h>
#include "../inst/include/utils.hpp"

//' @title C++ infection process (squire transmission model)
//' @description this is an internal function, you should use the R interface
//' for type checking, \code{\link{infection_process_cpp}}
//' @param parameters a list of parameters from \code{\link{get_parameters}}
//' @param states a \code{\link[individual]{CategoricalVariable}}
//' @param discrete_age a \code{\link[individual]{IntegerVariable}}
//' @param exposure a \code{\link[individual]{TargetedEvent}}
//' @param dt size of time step
// [[Rcpp::export]]
Rcpp::XPtr<process_t> infection_process_cpp_internal(
    Rcpp::List parameters,
    Rcpp::XPtr<CategoricalVariable> states,
    Rcpp::XPtr<IntegerVariable> discrete_age,
    Rcpp::XPtr<TargetedEvent> exposure,
    const double dt
) {

  // the states we need to pull
  std::vector<std::string> inf_states = {"IMild", "IAsymp", "ICase"};

  // vectors we can build once
  std::vector<double> beta(17, 0.);
  std::vector<double> lambda(17, 0.);

  return Rcpp::XPtr<process_t>(
    new process_t([parameters, states, discrete_age, exposure, dt, inf_states, beta, lambda](size_t t) mutable {

      // std::cout << " --- entering infection_process_cpp_internal on step: " << t << " --- \n";
      individual_index_t infectious = states->get_index_of(inf_states);
      // std::cout << "made infectious --- \n";

      if (infectious.size() > 0) {

        // current day (subtract one for zero-based indexing)
        size_t tnow = std::ceil((double)t * dt) - 1.;

        // group infection by age
        // std::cout << "group infection by age --- \n";
        std::vector<int> ages = discrete_age->get_values(infectious);
        std::vector<int> inf_ages = tab_bins(ages, 17);

        // calculate FoI for each age group
        // std::cout << "calculate FoI for each age group --- \n";
        Rcpp::NumericMatrix m = get_contact_matrix_cpp(parameters["mix_mat_set"], 0);
        std::fill(beta.begin(), beta.end(), get_vector_cpp(parameters["beta_set"], tnow));
        std::vector<double> m_inf_ages = matrix_vec_mult_cpp(m, inf_ages);
        std::transform(beta.begin(), beta.end(), m_inf_ages.begin(), lambda.begin(), std::multiplies<double>());

        // transition from S to E
        // std::cout << "transition from S to E --- \n";
        individual_index_t susceptible = states->get_index_of("S");
        std::vector<int> sus_ages = discrete_age->get_values(susceptible);

        // FoI for each susceptible person
        // std::cout << "FoI for each susceptible person --- \n";
        std::vector<double> lambda_sus(sus_ages.size());
        for (auto i = 0u; i < sus_ages.size(); ++i) {
          lambda_sus[i] = lambda[sus_ages[i] - 1] * dt;
          lambda_sus[i] = Rf_pexp(lambda_sus[i], 1., 1, 0);
        }

        // infected
        // std::cout << "infected --- \n";
        bitset_sample_multi_internal(susceptible, lambda_sus.begin(), lambda_sus.end());

        // newly infected queue the exposure event
        // std::cout << "newly infected queue the exposure event --- \n";
        if (susceptible.size() > 0) {
          // double zero{0.0};
          // exposure->schedule(susceptible, zero);
          exposure->schedule(susceptible, 0.0);
        }
      }

    }),
    true
  );
};
