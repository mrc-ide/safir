/* --------------------------------------------------------------------------------
 *  vaccination process for nimue vaccine model
 *  Sean L. Wu (slwood89@gmail.com)
 *  June 2021
 -------------------------------------------------------------------------------- */

#include <Rcpp.h>
#include <individual.h>
#include "../inst/include/utils.hpp"

//' @title C++ infection process (nimue vaccine model)
//' @description this is an internal function, you should use the R interface
//' for type checking, \code{\link{infection_process_cpp}}
//' @param parameters a list of parameters from \code{\link{get_parameters_nimue}}
//' @param states a \code{\link[individual]{CategoricalVariable}}
//' @param vaccine_states a \code{\link[individual]{IntegerVariable}}
//' @param discrete_age a \code{\link[individual]{IntegerVariable}}
//' @param exposure a \code{\link[individual]{TargetedEvent}}
//' @param dt size of time step
// [[Rcpp::export]]
Rcpp::XPtr<process_t> vaccination_process_nimue_cpp_internal(
    Rcpp::List parameters,
    Rcpp::XPtr<individual_index_t> eligible,
    Rcpp::XPtr<individual_index_t> vaccinated,
    Rcpp::XPtr<individual_index_t> empty,
    Rcpp::XPtr<IntegerVariable> discrete_age,
    Rcpp::XPtr<TargetedEvent> v0_to_v1v2,
    const double dt
) {

  // the states we need to pull
  std::vector<std::string> vaxx_states = {"S", "E", "R"};

  std::vector<double> pr(17);

  return Rcpp::XPtr<process_t>(
    new process_t([parameters, eligible, vaccinated, empty, discrete_age, v0_to_v1v2, dt, vaxx_states, pr](size_t t){

      // current day (subtract one for zero-based indexing)
      size_t tnow = std::ceil((double)t * dt) - 1.;

      double mv = get_vector_cpp(parameters["vaccine_set"], tnow);

      if (mv > 0.) {



      } // end if

    }),
    true
  );
};
