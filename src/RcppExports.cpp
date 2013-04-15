// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// comdistntInner
NumericMatrix comdistntInner(const int N, NumericMatrix dis, NumericMatrix x, bool exclude_conspecifics);
RcppExport SEXP sce_comdistntInner(SEXP NSEXP, SEXP disSEXP, SEXP xSEXP, SEXP exclude_conspecificsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    int N = Rcpp::as<int >(NSEXP);
    NumericMatrix dis = Rcpp::as<NumericMatrix >(disSEXP);
    NumericMatrix x = Rcpp::as<NumericMatrix >(xSEXP);
    bool exclude_conspecifics = Rcpp::as<bool >(exclude_conspecificsSEXP);
    NumericMatrix __result = comdistntInner(N, dis, x, exclude_conspecifics);
    return Rcpp::wrap(__result);
END_RCPP
}
// commdistInner
NumericMatrix commdistInner(const int N, NumericMatrix dis, NumericMatrix x);
RcppExport SEXP sce_commdistInner(SEXP NSEXP, SEXP disSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    int N = Rcpp::as<int >(NSEXP);
    NumericMatrix dis = Rcpp::as<NumericMatrix >(disSEXP);
    NumericMatrix x = Rcpp::as<NumericMatrix >(xSEXP);
    NumericMatrix __result = commdistInner(N, dis, x);
    return Rcpp::wrap(__result);
END_RCPP
}
