/* --------------------------------------------------------------------------------
 *  utilities for safir
 *  Sean L. Wu (slwood89@gmail.com)
 *  June 2021
 -------------------------------------------------------------------------------- */

#include "utils.hpp"

//' @title Cross tabulate two vectors with given margins
//' @description this is a replacement for \code{\link[base]{table}} that allows empty
//' cells because the margins have been specified. The input vectors \code{a} and \code{b}
//' must be the same length, this function does no argument checking.
//' @param a one set of observations
//' @param b another set of observations
//' @param a_margin margins for a (rows)
//' @param b_margin margins for b (cols)
//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix cross_tab_margins(
    const Rcpp::IntegerVector& a,
    const Rcpp::IntegerVector& b,
    const Rcpp::IntegerVector& a_margin,
    const Rcpp::IntegerVector& b_margin
) {

  auto r = a_margin.size();
  auto c = b_margin.size();
  Rcpp::IntegerMatrix out(r, c);

  for (auto i = 0u; i < a.size(); ++i) {
    out(a[i]-1, b[i]-1)++;
  }

  return out;
};
