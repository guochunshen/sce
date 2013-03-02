// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// nearest_neighbors
List nearest_neighbors(NumericVector xy, double rmax);
RcppExport SEXP sce_nearest_neighbors(SEXP xySEXP, SEXP rmaxSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    NumericVector xy = Rcpp::as<NumericVector >(xySEXP);
    double rmax = Rcpp::as<double >(rmaxSEXP);
    List __result = nearest_neighbors(xy, rmax);
    return Rcpp::wrap(__result);
END_RCPP
}
// rcpp_hello_world
List rcpp_hello_world();
RcppExport SEXP sce_rcpp_hello_world() {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    List __result = rcpp_hello_world();
    return Rcpp::wrap(__result);
END_RCPP
}
