/* --------------------------------------------------------------------------------
 *  utilities for safir (vaccine model, not nimue version)
 *  Sean L. Wu (slwood89@gmail.com)
 *  July 2021
 -------------------------------------------------------------------------------- */

#include "../inst/include/efficacy_vaccination.hpp"

//' @title Compute vaccine efficacy against infection from Ab titre
//' @param ab_titre a vector of Ab titres
//' @param parameters model parameters
//' @export
// [[Rcpp::export]]
std::vector<double> vaccine_efficacy_infection(
    const std::vector<double>& ab_titre,
    const Rcpp::List& parameters
) {
  // nt <- exp(ab_titre)
  std::vector<double> nt(ab_titre.size(), 0.0);
  std::transform(ab_titre.cbegin(), ab_titre.cend(), nt.begin(), [](const double nt_i) -> double {return std::exp(nt_i);});

  // ef_infection <- 1 / (1 + exp(-parameters$k * (log10(nt) - log10(parameters$ab_50)))) # reported efficacy in trials
  // ef_infection <- 1 - ef_infection (in the return line)
  double k = Rcpp::as<double>(parameters["k"]);
  double ab_50 = Rcpp::as<double>(parameters["ab_50"]);
  std::vector<double> ef_infection(ab_titre.size(), 0.0);
  std::transform(nt.cbegin(), nt.cend(), ef_infection.begin(), [k, ab_50](const double a) -> double {
    return 1.0 - (1.0 / (1.0 + std::exp(-k * (std::log10(a) - std::log10(ab_50)))));
  });

  // efficacy against infection
  return ef_infection;
};

//' @title Compute vaccine efficacy against severe disease from Ab titre
//' @description This needs the efficacy against infection because efficacy against severe disease,
//' conditional on breakthrough infection is what safir needs, which is computed as  1 - ((1 - efficacy_disease)/(1 - efficacy_infection)).
//' @param ab_titre a vector of Ab titres
//' @param ef_infection a vector of efficacy against infection from \code{\link{vaccine_efficacy_infection}}
//' @param parameters model parameters
//' @export
// [[Rcpp::export]]
std::vector<double> vaccine_efficacy_infection(
    const std::vector<double>& ab_titre,
    const std::vector<double>& ef_infection,
    const Rcpp::List& parameters
) {
  // nt <- exp(ab_titre)
  std::vector<double> nt(ab_titre.size(), 0.0);
  std::transform(ab_titre.cbegin(), ab_titre.cend(), nt.begin(), [](const double nt_i) -> double {return std::exp(nt_i);});

  // ef_severe_uncond <- 1 / (1 + exp(-parameters$k * (log10(nt) - log10(parameters$ab_50_severe))))
  // ef_severe <-  1 - ((1 - ef_severe_uncond)/(1 - ef_infection))
  // ef_severe <- 1 - ef_severe
  double k = Rcpp::as<double>(parameters["k"]);
  double ab_50_severe = Rcpp::as<double>(parameters["ab_50_severe"]);
  std::vector<double> ef_severe(ab_titre.size(), 0.0);
  std::transform(nt.cbegin(), nt.cend(), ef_infection.cbegin(), ef_severe.begin(), [k, ab_50_severe](const double nt_i, const double ef_inf_i) -> double {
    double ef_severe_uncond = 1.0 / (1.0 + std::exp(-k * (std::log10(nt_i) - std::log10(ab_50_severe))));
    double ef_severe = 1.0 - ((1.0 - ef_severe_uncond)/(1.0 - ef_inf_i));
    return 1.0 - ef_severe;
  });

  return ef_severe;
};
