/* --------------------------------------------------------------------------------
 *  utilities for safir
 *  Sean L. Wu (slwood89@gmail.com)
 *  June 2021
 -------------------------------------------------------------------------------- */

#ifndef UTILS_HPP
#define UTILS_HPP

#include <Rcpp.h>

Rcpp::IntegerMatrix cross_tab_margins(
    const Rcpp::IntegerVector& a,
    const Rcpp::IntegerVector& b,
    const int a_margin,
    const int b_margin
);

std::vector<int> tab_bins(
    const std::vector<int>& a,
    const int nbins
);

Rcpp::NumericMatrix get_mixing_matrix_cpp(
    SEXP array,
    const int i
);

#endif
