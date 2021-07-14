/* --------------------------------------------------------------------------------
 *  utilities for safir vaccination model (not nimue model)
 *  Sean L. Wu (slwood89@gmail.com)
 *  July 2021
 -------------------------------------------------------------------------------- */

#ifndef UTILS_VAXX_HPP
#define UTILS_VAXX_HPP

#include <Rcpp.h>
#include <individual.h>

std::vector<double> get_proportion_vaccinated_all_ages(
    const Rcpp::List variables,
    const int N_age,
    const int dose
);

#endif
