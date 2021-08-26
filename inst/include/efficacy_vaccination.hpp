/* --------------------------------------------------------------------------------
 *  efficacy functions for vaccination
 *  Sean L. Wu (slwood89@gmail.com)
 *  August 2021
 -------------------------------------------------------------------------------- */

#ifndef UTILS_EFFICACY_HPP
#define UTILS_EFFICACY_HPP

#include <Rcpp.h>
#include <individual.h>

//' @title Compute vaccine efficacy against infection from Ab titre
//' @param ab_titre a vector of Ab titres
//' @param parameters model parameters
std::vector<double> vaccine_efficacy_infection(
    const std::vector<double>& ab_titre,
    const Rcpp::List& parameters
);

//' @title Compute vaccine efficacy against severe disease from Ab titre
//' @description This needs the efficacy against infection because efficacy against severe disease,
//' conditional on breakthrough infection is what safir needs, which is computed as  1 - ((1 - efficacy_disease)/(1 - efficacy_infection)).
//' @param ab_titre a vector of Ab titres
//' @param ef_infection a vector of efficacy against infection from \code{\link{vaccine_efficacy_infection}}
//' @param parameters model parameters
std::vector<double> vaccine_efficacy_infection(
        const std::vector<double>& ab_titre,
        const std::vector<double>& ef_infection,
        const Rcpp::List& parameters
);

#endif
