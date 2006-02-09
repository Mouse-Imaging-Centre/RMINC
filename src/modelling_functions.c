#include "minc_reader.h"

SEXP t_test(SEXP voxel, SEXP grouping, SEXP n) {
  double *xvoxel, *xgrouping, *xn, *t;
  double x_mean, x_var, y_mean, y_var, x_sum, y_sum, x_n, y_n, s_p;
  int i;
  SEXP output;

  xgrouping = REAL(grouping);
  xvoxel = REAL(voxel);
  xn = REAL(n);

  PROTECT(output=allocVector(REALSXP, 1));
  t = REAL(output);

  /* compute sums and Ns */
  x_sum = 0;
  y_sum = 0;
  x_n = 0;
  y_n = 0;
  for (i=0; i < xn[0]; i++) {
    if (xgrouping[i] == 0) {
      x_n++;
      x_sum += xvoxel[i];
    }
    else if (xgrouping[i] == 1) {
      y_n++;
      y_sum += xvoxel[i];
    }
    else {
      error("Grouping value not 0 or 1\n");
    }
  }

  if (x_n == 0 || y_n == 0) {
    error("Each group must contain at least one subject\n");
  }

  x_mean = x_sum / x_n;
  y_mean = y_sum / y_n;

  x_var = 0;
  y_var = 0;

  /* compute the variances */
  for (i=0; i < xn[0]; i++) {
    if (xgrouping[i] == 0) {
      x_var += pow(xvoxel[i] - x_mean, 2);
    }
    else if (xgrouping[i] == 1) {
      y_var += pow(xvoxel[i] - y_mean, 2);
    }
  }
  x_var = x_var / x_n;
  y_var = y_var / y_n;

  /*
  Rprintf("Var (x,y) = %f %f\nMean (x,y) = %f %f\nN (x,y) %f %f\n", 
	  x_var, y_var, x_mean, y_mean, x_n, y_n);
  */

  //s_p = ( ((x_n - 1) * x_var) + ((y_n - 1) * y_var) ) / x_n + y_n - 2;
  t[0] = ( x_mean - y_mean ) / sqrt( (x_var / (x_n-1) ) + (y_var / (y_n-1)));

  UNPROTECT(1);
  return(output);
}


/* wilcoxon_rank_test: wilcoxon rank test (Mann-Whitney U test)
 * voxel: 1D array of doubles
 * grouping: 1D array of doubles
 * the test is against grouping == 0
 * NOTES:
 * does not handle ties in rankings correctly.
 */

SEXP wilcoxon_rank_test(SEXP voxel, SEXP grouping, SEXP na, SEXP nb) {
  double *xvoxel, *voxel_copy, *xgrouping, rank_sum, expected_rank_sum, *xW;
  double *xna, *xnb;
  int    *xindices;
  int *index;
  unsigned long    n, i;
  SEXP   indices, output;

  xgrouping = REAL(grouping);

  PROTECT(output=allocVector(REALSXP, 1));
  xW = REAL(output);
  xna = REAL(na);
  xnb = REAL(nb);

  //Rprintf("NA: %f NB: %f\n", *xna, *xnb);
  expected_rank_sum = *xna * *xnb + ((*xna * (*xna + 1)) / 2);

  n = LENGTH(voxel);
  xvoxel = REAL(voxel);
  voxel_copy = malloc(n * sizeof(double));
  index = malloc(n * sizeof(int));
  for (i=0; i < n; i++) {
    voxel_copy[i] = xvoxel[i];
    index[i] = i;
  }

  PROTECT(indices=allocVector(INTSXP, n));
  xindices = INTEGER(indices);

  rsort_with_index(voxel_copy, (int*) index, n);
  rank_sum = 0;
  
  for (i=0; i < n; i++) {
    //Rprintf("Index %d: %d %f\n", i, index[i]+1, xgrouping[i]);
    xindices[index[i]] = i+1;
    //if (xgrouping[i] == 0) 
    //expected_rank_sum += i+1;
  }
  for (i=0;  i< n; i++) {
    if (xgrouping[i] == 0)
      rank_sum += xindices[i];
  }
  //Rprintf("RANK SUM: %f\nEXPECTED SUM: %f\nW: %f\n", 
  //  rank_sum, expected_rank_sum, expected_rank_sum - rank_sum);
  xW[0] = expected_rank_sum - rank_sum;
  free(voxel_copy);
  free(index);
  UNPROTECT(2);
  return(output);
}
