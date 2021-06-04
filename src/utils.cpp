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
//' @param a_margin number of distinct values of a (rows)
//' @param b_margin number of distinct values of b (cols)
//' @examples
//' a <- 1:5
//' b <- c(1,2,3,1,2)
//' cross_tab_margins(a,b,5,3)
//' table(a,b)
//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix cross_tab_margins(
    const Rcpp::IntegerVector& a,
    const Rcpp::IntegerVector& b,
    const int a_margin,
    const int b_margin
) {

  Rcpp::IntegerMatrix out(a_margin, b_margin);

  for (auto i = 0u; i < a.size(); ++i) {
    out(a[i]-1, b[i]-1)++;
  }

  return out;
};
