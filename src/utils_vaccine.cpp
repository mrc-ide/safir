/* --------------------------------------------------------------------------------
 *  utilities for safir (vaccine model, not nimue version)
 *  Sean L. Wu (slwood89@gmail.com)
 *  July 2021
 -------------------------------------------------------------------------------- */

#include "../inst/include/utils_vaccine.hpp"

//' @title Get proportion of an age group that is vaccinated (multi-dose, no types)
//' @description Get proportion all age groups that have received a particular vaccine dose
//' by this timestep. This is similar to the function \code{\link[nimue]{coverage}} in the nimue package.
//' @param variables a list
//' @param N_age number of age groups
//' @param dose which dose to get proportion of age group who has received it (NOT 0-indexed)
//' @export
// [[Rcpp::export]]
std::vector<double> get_proportion_vaccinated_all_ages_cpp(
    const Rcpp::List variables,
    const int N_age,
    const int dose
){
  Rcpp::Environment discrete_age_env = Rcpp::as<Rcpp::Environment>(variables["discrete_age"]);
  Rcpp::XPtr<IntegerVariable> discrete_age = Rcpp::as<Rcpp::XPtr<IntegerVariable>>(discrete_age_env[".variable"]);
  Rcpp::List dose_time = Rcpp::as<Rcpp::List>(variables["dose_time"]);
  std::vector<double> out(N_age, 0.);
  for (int a = 1; a <= N_age; ++a) {
    individual_index_t age_bset = discrete_age->get_index_of_set(a);
    double N = age_bset.size();
    Rcpp::Environment dose_time_env = Rcpp::as<Rcpp::Environment>(dose_time[dose-1]);
    Rcpp::XPtr<IntegerVariable> dose_time_var = Rcpp::as<Rcpp::XPtr<IntegerVariable>>(dose_time_env[".variable"]);
    individual_index_t vaccinated_bset = dose_time_var->get_index_of_set(-1); // haven't gotten this dose
    vaccinated_bset = ~vaccinated_bset; // complement = vaccinated people
    vaccinated_bset &= age_bset;
    out[a-1] = static_cast<double>(vaccinated_bset.size()) / N;
  }
  return out;
};
