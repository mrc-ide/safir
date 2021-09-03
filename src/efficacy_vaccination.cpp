/* --------------------------------------------------------------------------------
 *  utilities for safir (vaccine model, not nimue version)
 *  Sean L. Wu (slwood89@gmail.com)
 *  July 2021
 -------------------------------------------------------------------------------- */

#include "../inst/include/efficacy_vaccination.hpp"

//' @title Compute vaccine efficacy against infection from Ab titre (C++)
//' @param ab_titre a vector of Ab titres
//' @param parameters model parameters
//' @return a numeric vector, 0 is maximally proective, 1 is maximally unprotective
//' @export
// [[Rcpp::export]]
std::vector<double> vaccine_efficacy_infection_cpp(
    const std::vector<double>& ab_titre,
    const Rcpp::List& parameters
) {
  double k = Rcpp::as<double>(parameters["k"]);
  double ab_50 = Rcpp::as<double>(parameters["ab_50"]);

  std::vector<double> ef_infection(ab_titre.size(), 1.0);

  for (auto i = 0u; i < ab_titre.size(); ++i) {
    if (std::isfinite(ab_titre[i])) {
      double nt = std::exp(ab_titre[i]);
      ef_infection[i] = 1.0 - (1.0 / (1.0 + std::exp(-k * (std::log10(nt) - std::log10(ab_50)))));
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
//' @param ab_titre a vector of Ab titres
//' @param ef_infection a vector of efficacy against infection from \code{\link{vaccine_efficacy_infection}}
//' @param parameters model parameters
//' @return a numeric vector, 0 is maximally proective, 1 is maximally unprotective
//' @export
// [[Rcpp::export]]
std::vector<double> vaccine_efficacy_severe_cpp(
    const std::vector<double>& ab_titre,
    const std::vector<double>& ef_infection,
    const Rcpp::List& parameters
) {
  double k = Rcpp::as<double>(parameters["k"]);
  double ab_50_severe = Rcpp::as<double>(parameters["ab_50_severe"]);

  std::vector<double> ef_severe(ab_titre.size(), 1.0);

  for (auto i = 0u; i < ab_titre.size(); ++i) {
    if (std::isfinite(ab_titre[i])) {
      double nt = std::exp(ab_titre[i]);
      double ef_severe_uncond = 1.0 / (1.0 + std::exp(-k * (std::log10(nt) - std::log10(ab_50_severe))));
      ef_severe[i] = 1.0 - ((1.0 - ef_severe_uncond) / ef_infection[i]); // 1 - (1 - ef_infection) goes from hazard reduction scale to efficacy, simplifies to ef_infection
      ef_severe[i] = 1.0 - ef_severe[i];
    } else {
      continue;
    }
  }

  // efficacy against severe disease
  return ef_severe;
};
