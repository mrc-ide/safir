/* --------------------------------------------------------------------------------
 *  utilities for safir (vaccine model, not nimue version)
 *  Sean L. Wu (slwood89@gmail.com)
 *  July 2021
 -------------------------------------------------------------------------------- */

#include "../inst/include/efficacy_nat.hpp"

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
