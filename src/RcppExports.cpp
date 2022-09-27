// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// pdDepthDists
double pdDepthDists(NumericVector new_dists, NumericMatrix dists, int N1);
RcppExport SEXP _depthtestr_pdDepthDists(SEXP new_distsSEXP, SEXP distsSEXP, SEXP N1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type new_dists(new_distsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type dists(distsSEXP);
    Rcpp::traits::input_parameter< int >::type N1(N1SEXP);
    rcpp_result_gen = Rcpp::wrap(pdDepthDists(new_dists, dists, N1));
    return rcpp_result_gen;
END_RCPP
}
// lcdDepthDists
double lcdDepthDists(NumericVector new_dists, NumericMatrix dists, int N1);
RcppExport SEXP _depthtestr_lcdDepthDists(SEXP new_distsSEXP, SEXP distsSEXP, SEXP N1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type new_dists(new_distsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type dists(distsSEXP);
    Rcpp::traits::input_parameter< int >::type N1(N1SEXP);
    rcpp_result_gen = Rcpp::wrap(lcdDepthDists(new_dists, dists, N1));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_depthtestr_pdDepthDists", (DL_FUNC) &_depthtestr_pdDepthDists, 3},
    {"_depthtestr_lcdDepthDists", (DL_FUNC) &_depthtestr_lcdDepthDists, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_depthtestr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
