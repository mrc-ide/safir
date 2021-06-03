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
    const Rcpp::IntegerVector& a_margin,
    const Rcpp::IntegerVector& b_margin
);

#endif
