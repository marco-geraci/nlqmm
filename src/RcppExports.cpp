// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// C_hfun_nlqmm
double C_hfun_nlqmm(NumericVector res, NumericVector u, NumericMatrix H, int N, int Q, double tau, double omicron);
RcppExport SEXP _nlqmm_C_hfun_nlqmm(SEXP resSEXP, SEXP uSEXP, SEXP HSEXP, SEXP NSEXP, SEXP QSEXP, SEXP tauSEXP, SEXP omicronSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type res(resSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type u(uSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type H(HSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type Q(QSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< double >::type omicron(omicronSEXP);
    rcpp_result_gen = Rcpp::wrap(C_hfun_nlqmm(res, u, H, N, Q, tau, omicron));
    return rcpp_result_gen;
END_RCPP
}
// C_loss_nlqmm
List C_loss_nlqmm(NumericVector res, NumericMatrix u, NumericVector weights, NumericMatrix Psiinv, double detH, double detPsi, int M, int N, int Q, double sigma, double tau, double omicron);
RcppExport SEXP _nlqmm_C_loss_nlqmm(SEXP resSEXP, SEXP uSEXP, SEXP weightsSEXP, SEXP PsiinvSEXP, SEXP detHSEXP, SEXP detPsiSEXP, SEXP MSEXP, SEXP NSEXP, SEXP QSEXP, SEXP sigmaSEXP, SEXP tauSEXP, SEXP omicronSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type res(resSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type u(uSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Psiinv(PsiinvSEXP);
    Rcpp::traits::input_parameter< double >::type detH(detHSEXP);
    Rcpp::traits::input_parameter< double >::type detPsi(detPsiSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type Q(QSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< double >::type omicron(omicronSEXP);
    rcpp_result_gen = Rcpp::wrap(C_loss_nlqmm(res, u, weights, Psiinv, detH, detPsi, M, N, Q, sigma, tau, omicron));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_nlqmm_C_hfun_nlqmm", (DL_FUNC) &_nlqmm_C_hfun_nlqmm, 7},
    {"_nlqmm_C_loss_nlqmm", (DL_FUNC) &_nlqmm_C_loss_nlqmm, 12},
    {NULL, NULL, 0}
};

RcppExport void R_init_nlqmm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
