#include <Rcpp.h>

#include "leapfrog-raw.h"

//' Simulate leapfrog model
//'
//' @param demp list of demographic input parameters (TODO: document)
//'
//'
//' @return a list containing projection output arrays
//'
//'  * `totpop1`: Total population by single-age, sex, and year.
//'  * `natdeaths`: Non-AIDS deaths to total population by single-age, sex, and year.
//' 
//' @details
//' The first year of `sx`, `asfr`, `srb`, and `netmig` is not used. This is assumed
//' to apply to the base year population (consistent with Spectrum).
//'
//' @export
//' 
// [[Rcpp::export]]
Rcpp::List
leapfrogR(const Rcpp::List& demp) {

  using namespace Rcpp;
  
  NumericVector Sx = demp["Sx"];
  Dimension d = Sx.attr("dim");
  const size_t proj_years = d[2];
  const int NG = 2;
  const int pAG = 81;
  const int pIDX_FERT = 15;
  const int pAG_FERT = 35;

  // allocate memory for return object
  NumericVector totpop1(pAG * NG * proj_years);
  totpop1.attr("dim") = NumericVector::create(pAG, NG, proj_years);

  NumericVector births(proj_years);

  NumericVector natdeaths(pAG * NG * proj_years);
  natdeaths.attr("dim") = NumericVector::create(pAG, NG, proj_years);

  
    leapfrog_sim<double, NG, pAG, pIDX_FERT, pAG_FERT>
      (REAL(demp["basepop"]),
       REAL(demp["Sx"]),
       REAL(demp["netmigr_adj"]),
       REAL(demp["asfr"]),
       REAL(demp["births_sex_prop"]),
       proj_years,
       REAL(totpop1),
       REAL(births),
       REAL(natdeaths));

  List ret = List::create(_("totpop1") = totpop1,
			  _("births") = births,
			  _("natdeaths") = natdeaths);
				      
  return ret;
}
