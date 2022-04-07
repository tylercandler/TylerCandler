#include <lapacke.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int main() {

  double A[9] = {4, -1, 0, -1, 4, -1, 0, -1, 4};
  double b[12] = {4, 5, 8, 2, 6, 1, 4, 7, 2, 6, 3, 7};

  lapack_int ipiv[3];
  int i;

  lapack_int n = 3;
  lapack_int nrhs = 4;
  lapack_int lda = 3;
  lapack_int ldb = 4;





lapack_int info = LAPACKE_dgesv(LAPACK_ROW_MAJOR,n,nrhs,A,lda,ipiv,b,ldb);

for (i=0;i<12;i++){
  printf("%lf\n",b[i]);
}


  return 0;
}






// -----------------------------------------------------------------------------
