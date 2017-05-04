#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double Rcpp_median(NumericVector x) {
  int lesser, greater, found = 0;
  int n = x.size();
  double med;
  for ( int i = 0; i < n; i++ ) {
    lesser = 0;
    greater = 0;
    for ( int j = 0; j < n; j++ ) {
      if ( x[j] < x[i] )
        lesser++;
      else if ( x[j] > x[i] )
        greater++;
    }
    if ( lesser == greater || abs(lesser - greater) == 1 ) {
      if ( n % 2 != 0 ) {
        med = x[i];
        break;
      }
      else {
        if ( found == 0 ) {
          found++;
          med = x[i];
        } else {
          med += x[i];
          med /= 2.0;
          break;
        }
      }
    }
  }
  return med;
}