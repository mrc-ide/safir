/* --------------------------------------------------------------------------------
 *  utilities for safir
 *  Sean L. Wu (slwood89@gmail.com)
 *  June 2021
 -------------------------------------------------------------------------------- */

#include "../inst/include/utils.hpp"

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

//' @title Tabulate a vector of observations
//' @description Tabulate a vector \code{a} whose values fall into a set of integers
//' of maximum value \code{nbins}. This function does no argument checking so please
//' ensure the maximum value of observations is not greater than \code{nbins}.
//' @param a a set of observations
//' @param nbins number of bins
//' @examples
//' nbin <- 10
//' a <- sample.int(n = nbin,size = 100,replace = T)
//' tabulate(bin = a,nbins = nbin)
//' tab_bins(a = a,nbins = nbin)
//' @export
// [[Rcpp::export]]
std::vector<int> tab_bins(
    const std::vector<int>& a,
    const int nbins
) {
  std::vector<int> out(nbins, 0);

  for (auto i = 0u; i < a.size(); ++i) {
    out[a[i]-1] += 1;
  }

  return out;
};

//' @title Get contact matrix
//' @description Get the contact matrix at some specific day (1st dimension of array).
//' @param array the mixing matrix array (days x age x age)
//' @param i the day (indexes the first dimension, assumes zero indexing)
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix get_contact_matrix_cpp(
    SEXP array,
    const int i
) {

  SEXP dims = Rf_getAttrib(array, R_DimSymbol);
  int d1 = INTEGER(dims)[0];
  int d2 = INTEGER(dims)[1];
  int d3 = INTEGER(dims)[2];

  if (i >= d1) {
    Rcpp::stop("bad value for i index passed to get_contact_matrix_cpp\n");
  }

  Rcpp::NumericMatrix out(d2,d3);

  double* array_ptr = REAL(array);

  for (auto j = 0u; j < d2; ++j) {
    for (auto k = 0u; k < d3; ++k) {
      out(j,k) = array_ptr[i + (j * d1) + (k * d1 * d2)];
    }
  }

  return out;
};

//' @title Get a beta value
//' @description Get a beta value at some specific day
//' @param beta_set the set of beta values
//' @param i the day (assumes zero indexing)
//' @export
// [[Rcpp::export]]
double get_beta_cpp(
    SEXP beta_set,
    const size_t i
){
  double* beta_set_ptr = REAL(beta_set);
  return beta_set_ptr[i];
};

//' @title Multiply a matrix by a integer vector
//' @param m a matrix
//' @param a a vector (must have length equal to number of columns of \code{m})
//' @export
// [[Rcpp::export]]
std::vector<double> matrix_vec_mult_cpp(
    const Rcpp::NumericMatrix& m,
    const std::vector<int>& a
) {
  std::vector<double> out(a.size(), 0.);

  for (auto j = 0u; j < m.ncol(); ++j) {
    for (auto i = 0u; i < m.ncol(); ++i) {
      out[i] += m(i,j) * (double)a[j];
    }
  }

  return out;
};

//' @title Multiply a matrix by a integer vector and a double vector
//' @param m a matrix
//' @param a a vector of int (must have length equal to number of columns of \code{m})
//' @param b a vector of double (must have length equal to number of columns of \code{m})
//' @export
// [[Rcpp::export]]
std::vector<double> matrix_2vec_mult_cpp(
    const Rcpp::NumericMatrix& m,
    const std::vector<int>& a,
    const std::vector<double>& b
) {
  std::vector<double> out(a.size(), 0.);

  for (auto j = 0u; j < m.ncol(); ++j) {
    for (auto i = 0u; i < m.ncol(); ++i) {
      out[i] += m(i,j) * (double)a[j] * b[j];
    }
  }

  return out;
};
