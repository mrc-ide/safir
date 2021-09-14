/* --------------------------------------------------------------------------------
 *  infection process for vaccination model (multiple doses, no types)
 *  Sean L. Wu (slwood89@gmail.com)
 *  July 2021
 -------------------------------------------------------------------------------- */

#include <Rcpp.h>
#include <individual.h>
#include "../inst/include/utils.hpp"
#include "../inst/include/efficacy_vaccination.hpp"

//' @title C++ infection process for vaccine model (multi-dose, no types)
//' @description this is an internal function, you should use the R interface
//' for type checking, \code{\link{infection_process_cpp}}
//' @param parameters a list of parameters from \code{\link{get_parameters}}
//' @param states a \code{\link[individual]{CategoricalVariable}}
//' @param discrete_age a \code{\link[individual]{IntegerVariable}}
//' @param ab_titre a \code{\link[individual]{DoubleVariable}}
//' @param exposure a \code{\link[individual]{TargetedEvent}}
//' @param dt size of time step
// [[Rcpp::export]]
Rcpp::XPtr<process_t> infection_process_vaccine_cpp_internal(
    Rcpp::List parameters,
    Rcpp::XPtr<CategoricalVariable> states,
    Rcpp::XPtr<IntegerVariable> discrete_age,
    Rcpp::XPtr<DoubleVariable> ab_titre,
    Rcpp::XPtr<TargetedEvent> exposure,
    const double dt
) {

  // the states we need to pull
  std::vector<std::string> inf_states = {"IMild", "IAsymp", "ICase"};

  // vectors we can build once
  std::vector<double> beta(17, 0.0);
  std::vector<double> lambda(17, 0.0);

  // parameters
  SEXP mix_mat_set = parameters["mix_mat_set"];
  SEXP beta_set = parameters["beta_set"];

  // infection process fn
  return Rcpp::XPtr<process_t>(
    new process_t([parameters, states, discrete_age, ab_titre, exposure, dt, inf_states, beta, lambda, mix_mat_set, beta_set](size_t t) mutable {

      individual_index_t infectious = states->get_index_of(inf_states);

      if (infectious.size() > 0) {

        // current day (subtract one for zero-based indexing)
        size_t tnow = std::ceil((double)t * dt) - 1.;

        // group infection by age
        std::vector<int> ages = discrete_age->get_values(infectious);
        std::vector<int> inf_ages = tab_bins(ages, 17);

        // calculate FoI for each age group
        Rcpp::NumericMatrix m = get_contact_matrix_cpp(mix_mat_set, 0);
        std::fill(beta.begin(), beta.end(), get_vector_cpp(beta_set, tnow));
        std::vector<double> m_inf_ages = matrix_vec_mult_cpp(m, inf_ages);
        std::transform(beta.begin(), beta.end(), m_inf_ages.begin(), lambda.begin(), std::multiplies<double>());

        // transition from S to E
        individual_index_t susceptible = states->get_index_of("S");
        std::vector<int> sus_ages = discrete_age->get_values(susceptible);

        // get efficacy
        std::vector<double> ab_titre_susceptible = ab_titre->get_values(susceptible);
        std::vector<double> infection_efficacy = vaccine_efficacy_infection_cpp(ab_titre_susceptible, parameters);

        // FoI for each susceptible person
        std::vector<double> lambda_sus(sus_ages.size());
        for (auto i = 0u; i < sus_ages.size(); ++i) {
          lambda_sus[i] = lambda[sus_ages[i] - 1] * infection_efficacy[i] * dt;
          lambda_sus[i] = Rf_pexp(lambda_sus[i], 1.0, 1, 0);
        }

        // infected
        bitset_sample_multi_internal(susceptible, lambda_sus.begin(), lambda_sus.end());

        // newly infected queue the exposure event
        if (susceptible.size() > 0) {
          exposure->schedule(susceptible, 0.0);
        }
      }

    }),
    true
  );
};
