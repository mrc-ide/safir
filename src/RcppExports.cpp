// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/safir_types.hpp"
#include <Rcpp.h>

using namespace Rcpp;

// infection_process_cpp_internal
Rcpp::XPtr<process_t> infection_process_cpp_internal(Rcpp::List parameters, Rcpp::XPtr<CategoricalVariable> states, Rcpp::XPtr<IntegerVariable> discrete_age, Rcpp::XPtr<TargetedEvent> exposure, const double dt);
RcppExport SEXP _safir_infection_process_cpp_internal(SEXP parametersSEXP, SEXP statesSEXP, SEXP discrete_ageSEXP, SEXP exposureSEXP, SEXP dtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type parameters(parametersSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<CategoricalVariable> >::type states(statesSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<IntegerVariable> >::type discrete_age(discrete_ageSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<TargetedEvent> >::type exposure(exposureSEXP);
    Rcpp::traits::input_parameter< const double >::type dt(dtSEXP);
    rcpp_result_gen = Rcpp::wrap(infection_process_cpp_internal(parameters, states, discrete_age, exposure, dt));
    return rcpp_result_gen;
END_RCPP
}
// cross_tab_margins
Rcpp::IntegerMatrix cross_tab_margins(const Rcpp::IntegerVector& a, const Rcpp::IntegerVector& b, const int a_margin, const int b_margin);
RcppExport SEXP _safir_cross_tab_margins(SEXP aSEXP, SEXP bSEXP, SEXP a_marginSEXP, SEXP b_marginSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const int >::type a_margin(a_marginSEXP);
    Rcpp::traits::input_parameter< const int >::type b_margin(b_marginSEXP);
    rcpp_result_gen = Rcpp::wrap(cross_tab_margins(a, b, a_margin, b_margin));
    return rcpp_result_gen;
END_RCPP
}
// tab_bins
std::vector<int> tab_bins(const std::vector<int>& a, const int nbins);
RcppExport SEXP _safir_tab_bins(SEXP aSEXP, SEXP nbinsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<int>& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const int >::type nbins(nbinsSEXP);
    rcpp_result_gen = Rcpp::wrap(tab_bins(a, nbins));
    return rcpp_result_gen;
END_RCPP
}
// get_contact_matrix_cpp
Rcpp::NumericMatrix get_contact_matrix_cpp(SEXP array, const int i);
RcppExport SEXP _safir_get_contact_matrix_cpp(SEXP arraySEXP, SEXP iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type array(arraySEXP);
    Rcpp::traits::input_parameter< const int >::type i(iSEXP);
    rcpp_result_gen = Rcpp::wrap(get_contact_matrix_cpp(array, i));
    return rcpp_result_gen;
END_RCPP
}
// get_beta_cpp
double get_beta_cpp(SEXP beta_set, const size_t i);
RcppExport SEXP _safir_get_beta_cpp(SEXP beta_setSEXP, SEXP iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type beta_set(beta_setSEXP);
    Rcpp::traits::input_parameter< const size_t >::type i(iSEXP);
    rcpp_result_gen = Rcpp::wrap(get_beta_cpp(beta_set, i));
    return rcpp_result_gen;
END_RCPP
}
// matrix_vec_mult_cpp
std::vector<double> matrix_vec_mult_cpp(const Rcpp::NumericMatrix& m, const std::vector<int> a);
RcppExport SEXP _safir_matrix_vec_mult_cpp(SEXP mSEXP, SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type m(mSEXP);
    Rcpp::traits::input_parameter< const std::vector<int> >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(matrix_vec_mult_cpp(m, a));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_safir_infection_process_cpp_internal", (DL_FUNC) &_safir_infection_process_cpp_internal, 5},
    {"_safir_cross_tab_margins", (DL_FUNC) &_safir_cross_tab_margins, 4},
    {"_safir_tab_bins", (DL_FUNC) &_safir_tab_bins, 2},
    {"_safir_get_contact_matrix_cpp", (DL_FUNC) &_safir_get_contact_matrix_cpp, 2},
    {"_safir_get_beta_cpp", (DL_FUNC) &_safir_get_beta_cpp, 2},
    {"_safir_matrix_vec_mult_cpp", (DL_FUNC) &_safir_matrix_vec_mult_cpp, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_safir(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
