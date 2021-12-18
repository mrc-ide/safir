/* --------------------------------------------------------------------------------
 *  infection process for vaccination model (multiple doses, no types)
 *  Sean L. Wu (slwood89@gmail.com)
 *  July 2021
 -------------------------------------------------------------------------------- */

#include <Rcpp.h>
#include <individual.h>
#include "../inst/include/utils.hpp"
#include "../inst/include/efficacy_nat.hpp"

// types for lambda functions used for calculations below
using get_inf_ages_func = std::function<std::vector<double>(const individual_index_t&, Rcpp::List, const size_t)>;

using calculate_nat_func = std::function<std::vector<double>(const individual_index_t&)>;


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
    Rcpp::List variables,
    Rcpp::XPtr<TargetedEvent> exposure,
    const double dt
) {

  // the states we need to pull
  std::vector<std::string> inf_states = {"IMild", "IAsymp", "ICase"};

  // vectors we can build once
  int N_age = Rcpp::as<int>(parameters["N_age"]);
  std::vector<double> beta(N_age, 0.0);
  std::vector<double> lambda_age(N_age, 0.0);

  // parameters
  SEXP mix_mat_set = parameters["mix_mat_set"];
  SEXP beta_set = parameters["beta_set"];
  SEXP lambda_external_vector = parameters["lambda_external"];

  // variables
  Rcpp::Environment discrete_age_R6 = Rcpp::as<Rcpp::Environment>(variables["discrete_age"]);
  Rcpp::XPtr<IntegerVariable> discrete_age(Rcpp::as<SEXP>(discrete_age_R6[".variable"]));

  Rcpp::Environment states_R6 = Rcpp::as<Rcpp::Environment>(variables["states"]);
  Rcpp::XPtr<CategoricalVariable> states(Rcpp::as<SEXP>(states_R6[".variable"]));

  Rcpp::Environment ab_titre_R6 = Rcpp::as<Rcpp::Environment>(variables["ab_titre"]);
  Rcpp::XPtr<DoubleVariable> ab_titre(Rcpp::as<SEXP>(ab_titre_R6[".variable"]));

  // calculate NAT
  calculate_nat_func calculate_nat;

  if (variables.containsElementNamed("ab_titre_inf")) {
    Rcpp::XPtr<DoubleVariable> ab_titre_inf(Rcpp::as<SEXP>(variables["ab_titre_inf"]));
    calculate_nat = [ab_titre, ab_titre_inf](const individual_index_t& index) -> std::vector<double> {

      std::vector<double> nat_vaccine = ab_titre->get_values(index);
      std::vector<double> nat_infection = ab_titre_inf->get_values(index);

      std::vector<double> nat(index.size());

      for (auto i = 0u; i < index.size(); ++i) {
        nat[i] = std::exp(nat_vaccine[i]) + std::exp(nat_infection[i]);
        nat[i] = std::log(nat[i]);
      }

      return nat;
    };
  } else {
    calculate_nat = [ab_titre](const individual_index_t& index) -> std::vector<double> {
      return ab_titre->get_values(index);
    };
  }

  // infection ages (weighted by NAT effect or not)
  get_inf_ages_func get_inf_ages;

  if (Rcpp::as<bool>(parameters["nt_efficacy_transmission"])) {
    get_inf_ages = [discrete_age, calculate_nat](const individual_index_t& infectious_bset, Rcpp::List parameters, const size_t timestep) -> std::vector<double> {
      int N_age = Rcpp::as<int>(parameters["N_age"]);
      std::vector<int> ages = discrete_age->get_values(infectious_bset);
      std::vector<double> nat_values = calculate_nat(infectious_bset);
      std::vector<double> inf_wt = vaccine_efficacy_transmission_cpp(nat_values, parameters, timestep);
      std::vector<double> inf_ages = tab_bins_weighted(ages, inf_wt, N_age);
      return(inf_ages);
    };
  } else {
    get_inf_ages = [discrete_age](const individual_index_t& infectious_bset, Rcpp::List parameters, const size_t timestep) -> std::vector<double> {
      int N_age = Rcpp::as<int>(parameters["N_age"]);
      std::vector<int> ages = discrete_age->get_values(infectious_bset);
      std::vector<double> inf_ages = tab_bins(ages, N_age);
      return(inf_ages);
    };
  }

  // infection process fn
  return Rcpp::XPtr<process_t>(
    new process_t([parameters, states, discrete_age, exposure, dt, inf_states, beta, lambda_age, mix_mat_set, beta_set, lambda_external_vector, get_inf_ages, calculate_nat](size_t t) mutable {

      // current day (subtract one for zero-based indexing)
      size_t day = std::ceil((double)t * dt) - 1.;

      // FoI from contact outside the population
      double lambda_external = get_vector_cpp(lambda_external_vector, day);

      // infectious classes
      individual_index_t infectious = states->get_index_of(inf_states);

      // susceptible persons
      individual_index_t susceptible = states->get_index_of("S");

      if (susceptible.size() > 0) {

        // FoI for each susceptible from external contacts
        std::vector<double> lambda(susceptible.size(), lambda_external);

        // FoI contribution from transmission
        if (infectious.size() > 0) {

          // group infectious persons by age
          std::vector<double> inf_ages = get_inf_ages(infectious, parameters, day);

          // calculate FoI on each susceptible age group
          Rcpp::NumericMatrix m = get_contact_matrix_cpp(mix_mat_set, 0);
          std::fill(beta.begin(), beta.end(), get_vector_cpp(beta_set, day));
          std::vector<double> m_inf_ages = matrix_vec_mult_cpp(m, inf_ages);
          std::transform(beta.begin(), beta.end(), m_inf_ages.begin(), lambda_age.begin(), std::multiplies<double>());

          // group susceptible persons by age
          std::vector<int> sus_ages = discrete_age->get_values(susceptible);

          // get vaccine efficacy
          std::vector<double> ab_titre_susceptible = calculate_nat(susceptible);
          std::vector<double> infection_efficacy = vaccine_efficacy_infection_cpp(ab_titre_susceptible, parameters, day);

          // FoI on each susceptible person from infectives
          for (auto i = 0u; i < sus_ages.size(); ++i) {
            lambda[i] += lambda_age[sus_ages[i] - 1] * infection_efficacy[i];
          }

        }

        // sample infection events in susceptible population
        std::transform(lambda.begin(), lambda.end(), lambda.begin(), [dt](const double l) -> double {
          return Rf_pexp(l * dt, 1.0, 1, 0);
        });
        bitset_sample_multi_internal(susceptible, lambda.begin(), lambda.end());

        // queue the exposure event
        if (susceptible.size() > 0) {
          exposure->schedule(susceptible, 0.0);
        }

      } // end if S>0

    }),
    true
  );
};
