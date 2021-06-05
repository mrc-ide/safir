/* --------------------------------------------------------------------------------
 *  utilities for safir
 *  Sean L. Wu (slwood89@gmail.com)
 *  June 2021
 -------------------------------------------------------------------------------- */

#ifndef UTILS_HPP
#define UTILS_HPP

#include <Rcpp.h>
#include <individual.h>

Rcpp::IntegerMatrix cross_tab_margins(
    const Rcpp::IntegerVector& a,
    const Rcpp::IntegerVector& b,
    const int a_margin,
    const int b_margin
);

Rcpp::NumericMatrix cross_tab_margins_internal(
        const std::vector<int>& a,
        const std::vector<int>& b,
        const int a_margin,
        const int b_margin
);

std::vector<int> tab_bins(
    const std::vector<int>& a,
    const int nbins
);

Rcpp::NumericMatrix get_contact_matrix_cpp(
    SEXP array,
    const int i
);

double get_vector_cpp(
    SEXP vector_set,
    const size_t i
);

std::vector<double> matrix_vec_mult_cpp(
    const Rcpp::NumericMatrix& m,
    const std::vector<int>& a
);

std::vector<double> matrix_2vec_mult_cpp(
    const Rcpp::NumericMatrix& m,
    const std::vector<double>& a,
    const std::vector<double>& b
);

std::vector<double> mult_2matrix_rowsum(
    const Rcpp::NumericMatrix& a,
    const Rcpp::NumericMatrix& b
);

double get_proportion_vaccinated_nimue_internal(
    Rcpp::XPtr<IntegerVariable> discrete_age,
    Rcpp::XPtr<individual_index_t> vaccinated,
    const int age
);

#endif
