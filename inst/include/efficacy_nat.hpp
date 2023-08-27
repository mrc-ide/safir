/* --------------------------------------------------------------------------------
 *  efficacy functions for vaccination
 *  Sean L. Wu (slwood89@gmail.com)
 *  August 2021
 -------------------------------------------------------------------------------- */

#ifndef UTILS_EFFICACY_HPP
#define UTILS_EFFICACY_HPP

#include <limits>

#include <Rcpp.h>
#include <individual.h>

using calculate_nat_func = std::function<std::vector<double>(const individual_index_t&, const size_t)>;

calculate_nat_func make_calculate_nat(
    const Rcpp::List& variables,
    const Rcpp::List& parameters
);

// used for testing only, calculate NAT
std::vector<double> test_make_calculate_nat_cpp(
    const Rcpp::List& variables,
    const Rcpp::List& parameters,
    Rcpp::XPtr<individual_index_t> index,
    const size_t day
);

std::vector<double> vaccine_efficacy_infection_cpp(
    const std::vector<double>& nat,
    const Rcpp::List& parameters,
    const size_t day
);

std::vector<double> vaccine_efficacy_severe_cpp(
    const std::vector<double>& nat,
    const std::vector<double>& ef_infection,
    const Rcpp::List& parameters,
    const size_t day
);

std::vector<double> vaccine_efficacy_transmission_cpp(
    const std::vector<double>& nat,
    const Rcpp::List& parameters,
    const size_t day
);

#endif
