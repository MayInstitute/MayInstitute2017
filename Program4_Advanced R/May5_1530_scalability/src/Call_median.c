#include <R.h>
#include <Rdefines.h>
#include <stdlib.h>

SEXP Call_median(SEXP x) {
  int lesser, greater, found = 0;
  int n = LENGTH(x);
  SEXP med;
  PROTECT(med = NEW_NUMERIC(1));
  double * pmed = REAL(med);
  double * px = REAL(x);
  for ( int i = 0; i < n; i++ ) {
    lesser = 0;
    greater = 0;
    for ( int j = 0; j < n; j++ ) {
      if ( px[j] < px[i] )
        lesser++;
      else if ( px[j] > px[i] )
        greater++;
    }
    if ( lesser == greater || abs(lesser - greater) == 1 ) {
      if ( n % 2 != 0 ) {
        *pmed = px[i];
        break;
      }
      else {
        if ( found == 0 ) {
          found++;
          *pmed = px[i];
        } else {
          *pmed += px[i];
          *pmed /= 2.0;
          break;
        }
      }
    }
  }
  UNPROTECT(1);
  return med;
}