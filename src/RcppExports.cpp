// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// wmwtest
extern SEXP wmwtest(SEXP matrix, SEXP labels, SEXP rtype);
RcppExport SEXP _bigde_wmwtest(SEXP matrixSEXP, SEXP labelsSEXP, SEXP rtypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type matrix(matrixSEXP);
    Rcpp::traits::input_parameter< SEXP >::type labels(labelsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type rtype(rtypeSEXP);
    rcpp_result_gen = Rcpp::wrap(wmwtest(matrix, labels, rtype));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_bigde_wmwtest", (DL_FUNC) &_bigde_wmwtest, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_bigde(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
