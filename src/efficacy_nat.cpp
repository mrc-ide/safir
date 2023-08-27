/* --------------------------------------------------------------------------------
 *  utilities for safir (vaccine model, not nimue version)
 *  Sean L. Wu (slwood89@gmail.com)
 *  July 2021
 -------------------------------------------------------------------------------- */

#include "../inst/include/efficacy_nat.hpp"

const double eps = std::numeric_limits<double>::epsilon();

calculate_nat_func make_calculate_nat(
    const Rcpp::List& variables,
    const Rcpp::List& parameters
) {

  // variables
  Rcpp::Environment ab_titre_R6 = Rcpp::as<Rcpp::Environment>(variables["ab_titre"]);
  Rcpp::XPtr<DoubleVariable> ab_titre(Rcpp::as<SEXP>(ab_titre_R6[".variable"]));

  calculate_nat_func calculate_nat;

  if (!variables.containsElementNamed("ab_titre_inf")) {
    Rcpp::stop("currently safir only correctly simulates models with seperate tracking of vaccine and infection derived NAT values");
  }

  Rcpp::Environment ab_titre_inf_R6 = Rcpp::as<Rcpp::Environment>(variables["ab_titre_inf"]);
  Rcpp::XPtr<DoubleVariable> ab_titre_inf(Rcpp::as<SEXP>(ab_titre_inf_R6[".variable"]));

  if (parameters.containsElementNamed("vp_time") && parameters["vp_time"] != R_NilValue) {

    Rcpp::Environment dose_time_R6 = Rcpp::as<Rcpp::Environment>(variables["dose_time"]);
    Rcpp::XPtr<IntegerVariable> dose_time(Rcpp::as<SEXP>(dose_time_R6[".variable"]));

    Rcpp::Environment dose_num_R6 = Rcpp::as<Rcpp::Environment>(variables["dose_num"]);
    Rcpp::XPtr<IntegerVariable> dose_num(Rcpp::as<SEXP>(dose_num_R6[".variable"]));

    calculate_nat = [ab_titre, ab_titre_inf, parameters, dose_time, dose_num](const individual_index_t& index, const size_t day) -> std::vector<double> {

      double vfr{1.0};
      if (parameters.containsElementNamed("vfr")) {
        SEXP vfr_sexp = parameters["vfr"];
        vfr = REAL(vfr_sexp)[day];
      }

      std::vector<double> nat_vaccine = ab_titre->get_values(index);
      std::vector<double> nat_infection = ab_titre_inf->get_values(index);

      double dt = Rcpp::as<double>(parameters["dt"]);

      SEXP vp_on_sexp = parameters["vp_time"];
      int* vp_on_ptr = INTEGER(vp_on_sexp);

      std::vector<int> dose_times = dose_time->get_values(index);
      std::vector<int> dose_nums = dose_num->get_values(index);

      std::vector<double> nat(index.size());

      double my_nat_infection, my_nat_vaccine;
      int my_dose_time;

      for (auto i = 0u; i < index.size(); ++i) {
        my_nat_infection = std::exp(nat_infection[i]);
        my_nat_vaccine = std::exp(nat_vaccine[i]);

        // nat_infection should always be scaled by vfr regardless
        my_nat_infection = std::max(eps, my_nat_infection / vfr);

        // if this person was vaccinated
        if (dose_nums[i] != 0) {
          // find individuals who were vaccinated when variant proof vaccine was on
          my_dose_time = static_cast<int>(std::ceil(static_cast<double>(dose_times[i]) * dt) - 1.0);
          if (vp_on_ptr[my_dose_time] == 0) {
            // apply vfr to those that were vaccinated not during variant proof window
            my_nat_vaccine = std::max(eps, my_nat_vaccine / vfr);
          }
        }

        // and combine these for overall NAT
        nat[i] = my_nat_infection + my_nat_vaccine;
      }

      return nat;
    };

  } else {

    calculate_nat = [ab_titre, ab_titre_inf, parameters](const individual_index_t& index, const size_t day) -> std::vector<double> {

      double vfr{1.0};
      if (parameters.containsElementNamed("vfr")) {
        SEXP vfr_sexp = parameters["vfr"];
        vfr = REAL(vfr_sexp)[day];
      }

      std::vector<double> nat_vaccine = ab_titre->get_values(index);
      std::vector<double> nat_infection = ab_titre_inf->get_values(index);

      std::vector<double> nat(index.size());

      for (auto i = 0u; i < index.size(); ++i) {
        nat_vaccine[i] = std::exp(nat_vaccine[i]);
        nat_infection[i] = std::exp(nat_infection[i]);
        nat[i] = nat_vaccine[i] + nat_infection[i];
        nat[i] = std::max(eps, nat[i] / vfr);
      }
      return nat;
    };

  }

  return calculate_nat;
};


//' @title used for testing only, calculate NAT
//' @param variables model variables
//' @param parameters model parameters
//' @param index a bitset
//' @param day current day
//' @return a numeric vector
//' @export
// [[Rcpp::export]]
std::vector<double> test_make_calculate_nat_cpp(
    const Rcpp::List& variables,
    const Rcpp::List& parameters,
    Rcpp::XPtr<individual_index_t> index,
    const size_t day
) {
  calculate_nat_func calc_nat = make_calculate_nat(variables, parameters);
  return calc_nat(*(index), day);
};


//' @title Compute vaccine efficacy against infection from Ab titre (C++)
//' @param nat a vector of NAT values
//' @param parameters model parameters
//' @param day current day
//' @return a numeric vector, 0 is maximally proective, 1 is maximally unprotective
//' @export
// [[Rcpp::export]]
std::vector<double> vaccine_efficacy_infection_cpp(
    const std::vector<double>& nat,
    const Rcpp::List& parameters,
    const size_t day
) {
  double k = Rcpp::as<double>(parameters["k"]);
  double ab_50 = Rcpp::as<double>(parameters["ab_50"]);
  double log10_ab_50 = std::log10(ab_50);

  std::vector<double> ef_infection(nat.size(), 1.0);
  double nt;

  for (auto i = 0u; i < nat.size(); ++i) {
    nt = nat[i];
    if (nt > 0.0) {
      ef_infection[i] = 1.0 - (1.0 / (1.0 + std::exp(-k * (std::log10(nt) - log10_ab_50))));
    } else {
      continue;
    }
  }

  // efficacy against infection
  return ef_infection;
};

//' @title Compute vaccine efficacy against severe disease from Ab titre (C++)
//' @description This needs the efficacy against infection because efficacy against severe disease,
//' conditional on breakthrough infection is what safir needs, which is computed as  1 - ((1 - efficacy_disease)/(1 - efficacy_infection)).
//' @param nat a vector of NAT values
//' @param ef_infection a vector of efficacy against infection from \code{\link{vaccine_efficacy_infection}}
//' @param parameters model parameters
//' @param day current day
//' @return a numeric vector, 0 is maximally proective, 1 is maximally unprotective
//' @export
// [[Rcpp::export]]
std::vector<double> vaccine_efficacy_severe_cpp(
    const std::vector<double>& nat,
    const std::vector<double>& ef_infection,
    const Rcpp::List& parameters,
    const size_t day
) {
  double k = Rcpp::as<double>(parameters["k"]);
  double ab_50_severe = Rcpp::as<double>(parameters["ab_50_severe"]);
  double log10_ab_50_severe = std::log10(ab_50_severe);

  std::vector<double> ef_severe(nat.size(), 1.0);
  double nt;

  for (auto i = 0u; i < nat.size(); ++i) {
    nt = nat[i];
    if (nt > 0.0) {
      double ef_severe_uncond = 1.0 / (1.0 + std::exp(-k * (std::log10(nt) - log10_ab_50_severe)));
      ef_severe[i] = 1.0 - ((1.0 - ef_severe_uncond) / ef_infection[i]); // 1 - (1 - ef_infection) goes from hazard reduction scale to efficacy, simplifies to ef_infection
      ef_severe[i] = 1.0 - ef_severe[i];
    } else {
      continue;
    }
  }

  // efficacy against severe disease
  return ef_severe;
};

//' @title Compute vaccine efficacy against onward transmission from Ab titre (C++)
//' @param nat a vector of NAT values
//' @param parameters model parameters.
//' @param day current day
//' @return a numeric vector, 0 is maximally protective, 1 is maximally unprotective
//' @export
// [[Rcpp::export]]
std::vector<double> vaccine_efficacy_transmission_cpp(
    const std::vector<double>& nat,
    const Rcpp::List& parameters,
    const size_t day
) {
  double k = Rcpp::as<double>(parameters["k"]);
  double ab_50 = Rcpp::as<double>(parameters["ab_50"]);
  double log10_ab_50 = std::log10(ab_50);
  double nt_transmission_factor = Rcpp::as<double>(parameters["nt_transmission_factor"]);

  std::vector<double> ef_transmission(nat.size(), 1.0);
  double nt;

  for (auto i = 0u; i < nat.size(); ++i) {
    nt = nat[i];
    if (nt > 0.0) {
      ef_transmission[i] = 1.0 - (1.0 / (1.0 + std::exp(-k * (std::log10(nt / nt_transmission_factor) - log10_ab_50))));
    } else {
      continue;
    }
  }

  // efficacy against infection
  return ef_transmission;
};
