#include <Rcpp.h>
#include <RcppEigen.h>

#include "leapfrog.h"


//' Simulate leapfrog model
//'
//' @param basepop base population: matrix indexed by (age, sex)
//' @param sx three-dimensional array of survival probabilities indexed by (age, sex, year)
//' @param asfr two-dimensional array of age-specific fertility rates (age, year)
//'
//'
//' @details
//' The first year of `sx`, `asfr`, `srb`, and `netmig` is not used. This is assumed
//' to apply to the base year population (consistent with Spectrum).
//'
//' @export
//' 
// [[Rcpp::export]]
Rcpp::NumericVector
leapfrogR(const Rcpp::NumericMatrix& basepop,
	  const Rcpp::NumericVector& Sx,
	  const Rcpp::NumericVector& netmigr,
	  const Rcpp::NumericVector& asfr,
	  const Rcpp::NumericVector& births_sex_prop) {

  Rcpp::Dimension d = Sx.attr("dim");
  const size_t proj_years = d[2];
  const int NG = 2;
  const int pAG = 81;
  const int pIDX_FERT = 15;
  const int pAG_FERT = 35;

  // return object
  Rcpp::NumericVector pop(pAG * NG * proj_years);
  pop.attr("dim") = Rcpp::NumericVector::create(pAG, NG, proj_years);

  leapfrog_sim<double, NG, pAG, pIDX_FERT, pAG_FERT>
    (REAL(basepop),
     REAL(Sx),
     REAL(netmigr),
     REAL(asfr),
     REAL(births_sex_prop),
     proj_years,
     REAL(pop));
  
  return pop;
}


