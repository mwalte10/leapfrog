// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// make_leslie_matrixR
Eigen::SparseMatrix<double> make_leslie_matrixR(const Eigen::Map<Eigen::ArrayXd> sx, const Eigen::Map<Eigen::ArrayXd> fx, double srb, double age_span, int fx_idx);
RcppExport SEXP _leapfrog_make_leslie_matrixR(SEXP sxSEXP, SEXP fxSEXP, SEXP srbSEXP, SEXP age_spanSEXP, SEXP fx_idxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::ArrayXd> >::type sx(sxSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::ArrayXd> >::type fx(fxSEXP);
    Rcpp::traits::input_parameter< double >::type srb(srbSEXP);
    Rcpp::traits::input_parameter< double >::type age_span(age_spanSEXP);
    Rcpp::traits::input_parameter< int >::type fx_idx(fx_idxSEXP);
    rcpp_result_gen = Rcpp::wrap(make_leslie_matrixR(sx, fx, srb, age_span, fx_idx));
    return rcpp_result_gen;
END_RCPP
}
// ccmppR
Eigen::MatrixXd ccmppR(const Eigen::Map<Eigen::VectorXd> basepop, const Eigen::Map<Eigen::MatrixXd> sx, const Eigen::Map<Eigen::MatrixXd> fx, const Eigen::Map<Eigen::VectorXd> srb, double age_span, int fx_idx);
RcppExport SEXP _leapfrog_ccmppR(SEXP basepopSEXP, SEXP sxSEXP, SEXP fxSEXP, SEXP srbSEXP, SEXP age_spanSEXP, SEXP fx_idxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd> >::type basepop(basepopSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type sx(sxSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type fx(fxSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd> >::type srb(srbSEXP);
    Rcpp::traits::input_parameter< double >::type age_span(age_spanSEXP);
    Rcpp::traits::input_parameter< int >::type fx_idx(fx_idxSEXP);
    rcpp_result_gen = Rcpp::wrap(ccmppR(basepop, sx, fx, srb, age_span, fx_idx));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_leapfrog_make_leslie_matrixR", (DL_FUNC) &_leapfrog_make_leslie_matrixR, 5},
    {"_leapfrog_ccmppR", (DL_FUNC) &_leapfrog_ccmppR, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_leapfrog(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
