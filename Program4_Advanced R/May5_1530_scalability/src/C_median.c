#include <R.h>
#include <stdlib.h>

void C_median(double * x, int * n, double * med) {
  int lesser, greater, found = 0;
  for ( int i = 0; i < *n; i++ ) {
    lesser = 0;
    greater = 0;
    for ( int j = 0; j < *n; j++ ) {
      if ( x[j] < x[i] )
        lesser++;
      else if ( x[j] > x[i] )
        greater++;
    }
    if ( lesser == greater || abs(lesser - greater) == 1 ) {
      if ( *n % 2 != 0 ) {
        *med = x[i];
        break;
      }
      else {
        if ( found == 0 ) {
          found++;
          *med = x[i];
        } else {
          *med += x[i];
          *med /= 2.0;
          break;
        }
      }
    }
  }
}
