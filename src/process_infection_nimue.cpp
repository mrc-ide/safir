/* --------------------------------------------------------------------------------
 *  infection process for nimue vaccine model
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
//' @param vaccine_states a \code{\link[individual]{IntegerVariable}}
//' @param discrete_age a \code{\link[individual]{IntegerVariable}}
//' @param exposure a \code{\link[individual]{TargetedEvent}}
//' @param dt size of time step
// [[Rcpp::export]]
Rcpp::XPtr<process_t> infection_process_nimue_cpp_internal(
    Rcpp::List parameters,
    Rcpp::XPtr<CategoricalVariable> states,
    Rcpp::XPtr<IntegerVariable> vaccine_states,
    Rcpp::XPtr<IntegerVariable> discrete_age,
    Rcpp::XPtr<TargetedEvent> exposure,
    const double dt
) {

  // the states we need to pull
  std::vector<std::string> inf_states = {"IMild", "IAsymp", "ICase"};

  return Rcpp::XPtr<process_t>(
    new process_t([parameters, states, vaccine_states, discrete_age, exposure, dt, inf_states](size_t t){

      individual_index_t infectious = states->get_index_of(inf_states);

      if (infectious.size() > 0) {

        // current day (subtract one for zero-based indexing)
        size_t tnow = std::ceil((double)t * dt) - 1.;

        // infection by vaccine status
        std::vector<int> inf_vaxx = vaccine_states->get_values(infectious);

        // infection by age
        std::vector<int> ages = discrete_age->get_values(infectious);

        // compute cross tab for relative infectiousness, multiply by that matrix, and sum it out
        Rcpp::NumericMatrix inf_age_vax = cross_tab_margins_internal(ages, inf_vaxx, 17, 4);
        std::vector<double> inf_ages = mult_2matrix_rowsum(inf_age_vax, Rcpp::as<Rcpp::NumericMatrix>(parameters["rel_infectiousness_vaccinated"]));

        // calculate FoI for each age group
        Rcpp::NumericMatrix m = get_contact_matrix_cpp(parameters["mix_mat_set"], 0);
        std::vector<double> beta(17, get_vector_cpp(parameters["beta_set"], tnow));
        std::vector<double> m_inf_ages_rel = matrix_2vec_mult_cpp(m, inf_ages, Rcpp::as<std::vector<double>>(parameters["rel_infectiousness"]));
        std::vector<double> lambda(17);
        std::transform(beta.begin(), beta.end(), m_inf_ages_rel.begin(), lambda.begin(), std::multiplies<double>());

        // FoI for each susceptible person
        individual_index_t susceptible = states->get_index_of("S");
        std::vector<int> sus_vaxx = vaccine_states->get_values(susceptible);
        std::vector<int> sus_ages = discrete_age->get_values(susceptible);

        SEXP vaccine_efficacy_infection = parameters["vaccine_efficacy_infection"];
        SEXP dims = Rf_getAttrib(vaccine_efficacy_infection, R_DimSymbol);
        int d1 = INTEGER(dims)[0];
        int d2 = INTEGER(dims)[1];
        int d3 = INTEGER(dims)[2];
        double* vaccine_efficacy_infection_ptr = REAL(vaccine_efficacy_infection);

        std::vector<double> lambda_sus(sus_ages.size());

        for (auto i = 0u; i < sus_ages.size(); ++i) {
          int age = sus_ages[i] - 1;
          int vax = sus_vaxx[i] - 1;
          lambda_sus[i] = lambda[age] * vaccine_efficacy_infection_ptr[tnow + (age * d1) + (vax * d1 * d2)] * dt;
          lambda_sus[i] = Rf_pexp(lambda_sus[i], 1., 1, 0);
        }

        // infected
        bitset_sample_multi_internal(susceptible, lambda_sus.begin(), lambda_sus.end());

        // newly infected queue the exposure event
        if (susceptible.size() > 0) {
          exposure->schedule(susceptible, 0.);
        }
      }

    }),
    true
  );
};
