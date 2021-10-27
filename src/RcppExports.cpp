// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/safir_types.hpp"
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// vaccine_efficacy_infection_cpp
std::vector<double> vaccine_efficacy_infection_cpp(const std::vector<double>& ab_titre, const Rcpp::List& parameters);
RcppExport SEXP _safir_vaccine_efficacy_infection_cpp(SEXP ab_titreSEXP, SEXP parametersSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double>& >::type ab_titre(ab_titreSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type parameters(parametersSEXP);
    rcpp_result_gen = Rcpp::wrap(vaccine_efficacy_infection_cpp(ab_titre, parameters));
    return rcpp_result_gen;
END_RCPP
}
// vaccine_efficacy_severe_cpp
std::vector<double> vaccine_efficacy_severe_cpp(const std::vector<double>& ab_titre, const std::vector<double>& ef_infection, const Rcpp::List& parameters);
RcppExport SEXP _safir_vaccine_efficacy_severe_cpp(SEXP ab_titreSEXP, SEXP ef_infectionSEXP, SEXP parametersSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double>& >::type ab_titre(ab_titreSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type ef_infection(ef_infectionSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type parameters(parametersSEXP);
    rcpp_result_gen = Rcpp::wrap(vaccine_efficacy_severe_cpp(ab_titre, ef_infection, parameters));
    return rcpp_result_gen;
END_RCPP
}
// vaccination_process_nimue_cpp_internal
Rcpp::XPtr<process_t> vaccination_process_nimue_cpp_internal(Rcpp::List parameters, Rcpp::XPtr<CategoricalVariable> states, Rcpp::XPtr<individual_index_t> eligible, Rcpp::XPtr<individual_index_t> vaccinated, Rcpp::XPtr<individual_index_t> empty, Rcpp::XPtr<IntegerVariable> discrete_age, Rcpp::XPtr<TargetedEvent> v0_to_v1v2, const double dt);
RcppExport SEXP _safir_vaccination_process_nimue_cpp_internal(SEXP parametersSEXP, SEXP statesSEXP, SEXP eligibleSEXP, SEXP vaccinatedSEXP, SEXP emptySEXP, SEXP discrete_ageSEXP, SEXP v0_to_v1v2SEXP, SEXP dtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type parameters(parametersSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<CategoricalVariable> >::type states(statesSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<individual_index_t> >::type eligible(eligibleSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<individual_index_t> >::type vaccinated(vaccinatedSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<individual_index_t> >::type empty(emptySEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<IntegerVariable> >::type discrete_age(discrete_ageSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<TargetedEvent> >::type v0_to_v1v2(v0_to_v1v2SEXP);
    Rcpp::traits::input_parameter< const double >::type dt(dtSEXP);
    rcpp_result_gen = Rcpp::wrap(vaccination_process_nimue_cpp_internal(parameters, states, eligible, vaccinated, empty, discrete_age, v0_to_v1v2, dt));
    return rcpp_result_gen;
END_RCPP
}
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
// infection_process_nimue_cpp_internal
Rcpp::XPtr<process_t> infection_process_nimue_cpp_internal(Rcpp::List parameters, Rcpp::XPtr<CategoricalVariable> states, Rcpp::XPtr<IntegerVariable> vaccine_states, Rcpp::XPtr<IntegerVariable> discrete_age, Rcpp::XPtr<TargetedEvent> exposure, const double dt);
RcppExport SEXP _safir_infection_process_nimue_cpp_internal(SEXP parametersSEXP, SEXP statesSEXP, SEXP vaccine_statesSEXP, SEXP discrete_ageSEXP, SEXP exposureSEXP, SEXP dtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type parameters(parametersSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<CategoricalVariable> >::type states(statesSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<IntegerVariable> >::type vaccine_states(vaccine_statesSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<IntegerVariable> >::type discrete_age(discrete_ageSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<TargetedEvent> >::type exposure(exposureSEXP);
    Rcpp::traits::input_parameter< const double >::type dt(dtSEXP);
    rcpp_result_gen = Rcpp::wrap(infection_process_nimue_cpp_internal(parameters, states, vaccine_states, discrete_age, exposure, dt));
    return rcpp_result_gen;
END_RCPP
}
// infection_process_vaccine_cpp_internal
Rcpp::XPtr<process_t> infection_process_vaccine_cpp_internal(Rcpp::List parameters, Rcpp::XPtr<CategoricalVariable> states, Rcpp::XPtr<IntegerVariable> discrete_age, Rcpp::XPtr<DoubleVariable> ab_titre, Rcpp::XPtr<TargetedEvent> exposure, const double dt);
RcppExport SEXP _safir_infection_process_vaccine_cpp_internal(SEXP parametersSEXP, SEXP statesSEXP, SEXP discrete_ageSEXP, SEXP ab_titreSEXP, SEXP exposureSEXP, SEXP dtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type parameters(parametersSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<CategoricalVariable> >::type states(statesSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<IntegerVariable> >::type discrete_age(discrete_ageSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<DoubleVariable> >::type ab_titre(ab_titreSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<TargetedEvent> >::type exposure(exposureSEXP);
    Rcpp::traits::input_parameter< const double >::type dt(dtSEXP);
    rcpp_result_gen = Rcpp::wrap(infection_process_vaccine_cpp_internal(parameters, states, discrete_age, ab_titre, exposure, dt));
    return rcpp_result_gen;
END_RCPP
}
// compare_floats
bool compare_floats(const double a, const double b);
RcppExport SEXP _safir_compare_floats(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(compare_floats(a, b));
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
// cross_tab_doses_age
Rcpp::IntegerMatrix cross_tab_doses_age(Rcpp::XPtr<IntegerVariable> doses, Rcpp::XPtr<IntegerVariable> age, const size_t num_doses, const size_t num_ages);
RcppExport SEXP _safir_cross_tab_doses_age(SEXP dosesSEXP, SEXP ageSEXP, SEXP num_dosesSEXP, SEXP num_agesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<IntegerVariable> >::type doses(dosesSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<IntegerVariable> >::type age(ageSEXP);
    Rcpp::traits::input_parameter< const size_t >::type num_doses(num_dosesSEXP);
    Rcpp::traits::input_parameter< const size_t >::type num_ages(num_agesSEXP);
    rcpp_result_gen = Rcpp::wrap(cross_tab_doses_age(doses, age, num_doses, num_ages));
    return rcpp_result_gen;
END_RCPP
}
// cross_tab_margins_internal
Rcpp::NumericMatrix cross_tab_margins_internal(const std::vector<int>& a, const std::vector<int>& b, const int a_margin, const int b_margin);
RcppExport SEXP _safir_cross_tab_margins_internal(SEXP aSEXP, SEXP bSEXP, SEXP a_marginSEXP, SEXP b_marginSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<int>& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const std::vector<int>& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const int >::type a_margin(a_marginSEXP);
    Rcpp::traits::input_parameter< const int >::type b_margin(b_marginSEXP);
    rcpp_result_gen = Rcpp::wrap(cross_tab_margins_internal(a, b, a_margin, b_margin));
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
// get_vector_cpp
double get_vector_cpp(SEXP vector_set, const size_t i);
RcppExport SEXP _safir_get_vector_cpp(SEXP vector_setSEXP, SEXP iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type vector_set(vector_setSEXP);
    Rcpp::traits::input_parameter< const size_t >::type i(iSEXP);
    rcpp_result_gen = Rcpp::wrap(get_vector_cpp(vector_set, i));
    return rcpp_result_gen;
END_RCPP
}
// matrix_vec_mult_cpp
std::vector<double> matrix_vec_mult_cpp(const Rcpp::NumericMatrix& m, const std::vector<int>& a);
RcppExport SEXP _safir_matrix_vec_mult_cpp(SEXP mSEXP, SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type m(mSEXP);
    Rcpp::traits::input_parameter< const std::vector<int>& >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(matrix_vec_mult_cpp(m, a));
    return rcpp_result_gen;
END_RCPP
}
// matrix_2vec_mult_cpp
std::vector<double> matrix_2vec_mult_cpp(const Rcpp::NumericMatrix& m, const std::vector<double>& a, const std::vector<double>& b);
RcppExport SEXP _safir_matrix_2vec_mult_cpp(SEXP mSEXP, SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type m(mSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(matrix_2vec_mult_cpp(m, a, b));
    return rcpp_result_gen;
END_RCPP
}
// mult_2matrix_rowsum
std::vector<double> mult_2matrix_rowsum(const Rcpp::NumericMatrix& a, const Rcpp::NumericMatrix& b);
RcppExport SEXP _safir_mult_2matrix_rowsum(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(mult_2matrix_rowsum(a, b));
    return rcpp_result_gen;
END_RCPP
}
// get_proportion_vaccinated_nimue_internal
double get_proportion_vaccinated_nimue_internal(Rcpp::XPtr<IntegerVariable> discrete_age, Rcpp::XPtr<individual_index_t> vaccinated, const int age);
RcppExport SEXP _safir_get_proportion_vaccinated_nimue_internal(SEXP discrete_ageSEXP, SEXP vaccinatedSEXP, SEXP ageSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<IntegerVariable> >::type discrete_age(discrete_ageSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<individual_index_t> >::type vaccinated(vaccinatedSEXP);
    Rcpp::traits::input_parameter< const int >::type age(ageSEXP);
    rcpp_result_gen = Rcpp::wrap(get_proportion_vaccinated_nimue_internal(discrete_age, vaccinated, age));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_safir_vaccine_efficacy_infection_cpp", (DL_FUNC) &_safir_vaccine_efficacy_infection_cpp, 2},
    {"_safir_vaccine_efficacy_severe_cpp", (DL_FUNC) &_safir_vaccine_efficacy_severe_cpp, 3},
    {"_safir_vaccination_process_nimue_cpp_internal", (DL_FUNC) &_safir_vaccination_process_nimue_cpp_internal, 8},
    {"_safir_infection_process_cpp_internal", (DL_FUNC) &_safir_infection_process_cpp_internal, 5},
    {"_safir_infection_process_nimue_cpp_internal", (DL_FUNC) &_safir_infection_process_nimue_cpp_internal, 6},
    {"_safir_infection_process_vaccine_cpp_internal", (DL_FUNC) &_safir_infection_process_vaccine_cpp_internal, 6},
    {"_safir_compare_floats", (DL_FUNC) &_safir_compare_floats, 2},
    {"_safir_cross_tab_margins", (DL_FUNC) &_safir_cross_tab_margins, 4},
    {"_safir_cross_tab_doses_age", (DL_FUNC) &_safir_cross_tab_doses_age, 4},
    {"_safir_cross_tab_margins_internal", (DL_FUNC) &_safir_cross_tab_margins_internal, 4},
    {"_safir_tab_bins", (DL_FUNC) &_safir_tab_bins, 2},
    {"_safir_get_contact_matrix_cpp", (DL_FUNC) &_safir_get_contact_matrix_cpp, 2},
    {"_safir_get_vector_cpp", (DL_FUNC) &_safir_get_vector_cpp, 2},
    {"_safir_matrix_vec_mult_cpp", (DL_FUNC) &_safir_matrix_vec_mult_cpp, 2},
    {"_safir_matrix_2vec_mult_cpp", (DL_FUNC) &_safir_matrix_2vec_mult_cpp, 3},
    {"_safir_mult_2matrix_rowsum", (DL_FUNC) &_safir_mult_2matrix_rowsum, 2},
    {"_safir_get_proportion_vaccinated_nimue_internal", (DL_FUNC) &_safir_get_proportion_vaccinated_nimue_internal, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_safir(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
