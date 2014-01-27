#include "minc_reader.h"
#include "R_ext/Lapack.h"
#include "R_ext/Applic.h"

/* minimum computation to speed up bits of qvalue */

void qvalue_min(double *qvalue, int *u, int *length, double *out) {
  int i;
  
  //Rprintf("Length: %d\n", length[0]);

  length[0] -= 1;

  //Rprintf("Length: %d\n", length[0]);

  if (qvalue[u[length[0]]] > 1) {
    out[u[length[0]]-1] = 1;
  }
  else {
    out[u[length[0]]-1] = qvalue[u[length[0]]-1];
  }
  
  for (i = length[0]-1; i >= 0; i--) {
    //u[i] -= 1;
    //Rprintf("%f %f %d %d\n", qvalue[u[i]-1], qvalue[u[i+1]-1], i, u[i]-1);
    if ((qvalue[u[i]-1] > 1) && (qvalue[u[i+1]-1] > 1)) {
      out[u[i]-1] = 1;
      //Rprintf("out = %f\n", 1.0);
    }
    else if (qvalue[u[i]-1] < qvalue[u[i+1]-1]) {
      out[u[i]-1] = qvalue[u[i]-1];
      //Rprintf("out = %f\n", qvalue[u[i]-1]);
    }
    else {
      out[u[i]-1] = qvalue[u[i+1]-1];
      //Rprintf("out = %f\n", qvalue[u[i+1]-1]);
    }
  }
}

/* compute a paired t test given a voxel and grouping
 *
 * Key assumption: that the first voxel belonging to group 0 is to be
 * matched with the first voxel belonging to group 1, etc.
 */

SEXP paired_t_test(SEXP voxel, SEXP grouping) {
  double *xvoxel, *xgrouping, *t;
  double *group0, *group1;
  double mean_difference, sd;
  int i, n, n2, count0, count1;
  SEXP output;

  xvoxel = REAL(voxel);
  xgrouping = REAL(grouping);

  n = LENGTH(grouping);
  n2 = n/2;

  PROTECT(output=allocVector(REALSXP, 1));
  t = REAL(output);

  group0 = malloc(sizeof(double) * n2);
  group1 = malloc(sizeof(double) * n2);

  count0 = 0;
  count1 = 0;

  /* assign voxels to each groups vector */
  for (i=0; i < n; i++) {
    if (xgrouping[i] == 0) {
      group0[count0] = xvoxel[i];
      count0++;
    }
    else if (xgrouping[i] == 1) {
      group1[count1] = xvoxel[i];
      count1++;
    }
  }
  

  mean_difference = 0;
  for (i=0; i < n2; i++) {
    mean_difference += group0[i] - group1[i];
  }
  mean_difference = mean_difference / n2;
  
  sd = 0;
  for (i=0; i < n2; i++) {
    sd += pow((group0[i] - group1[i]) - mean_difference, 2);
  }
  sd = sqrt(sd/(n2-1));

  t[0] = mean_difference / (sd / sqrt(n2));

  free(group0);
  free(group1);
  UNPROTECT(1);
  return(output);
}
  

/* compute a  t test given a voxel and grouping */
SEXP t_test(SEXP voxel, SEXP grouping) {
  double *xvoxel, *xgrouping, *t;
  double x_mean, x_var, y_mean, y_var, x_sum, y_sum, x_n, y_n;
  int i, n;
  SEXP output;

  n = LENGTH(grouping);
  xgrouping = REAL(grouping);
  xvoxel = REAL(voxel);

  PROTECT(output=allocVector(REALSXP, 1));
  t = REAL(output);

  /* compute sums and Ns */
  x_sum = 0;
  y_sum = 0;
  x_n = 0;
  y_n = 0;
  for (i=0; i < n; i++) {
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
  for (i=0; i < n; i++) {
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

SEXP wilcoxon_rank_test(SEXP voxel, SEXP grouping) {
  double *xvoxel, *voxel_copy, *xgrouping, rank_sum, expected_rank_sum, *xW;
  int    *xindices;
  int *index;
  unsigned long    n, i, na, nb;
  SEXP   indices, output;

  n = LENGTH(grouping);
  xgrouping = REAL(grouping);

  PROTECT(output=allocVector(REALSXP, 1));
  xW = REAL(output);

  na = 0;
  nb = 0;
  for (i=0; i<n; i++) {
    if (xgrouping[i] == 0)
      na++;
    else if (xgrouping[i] == 1)
      nb++;
    else
      error("Each element in grouping must be either 0 or 1\n");
  }

  //Rprintf("NA: %f NB: %f\n", na, nb);
  expected_rank_sum = na * nb + ((na * (na + 1)) / 2);

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

SEXP voxel_sum(SEXP Svoxel, SEXP Sn_groups, SEXP Sgroupings) {
  double *voxel, *groupings, *n_groups;
  double *sum_voxel;
  int n_subjects, i;
  SEXP Soutput;


  voxel = REAL(Svoxel);
  n_groups = REAL(Sn_groups);
  groupings = REAL(Sgroupings);

  n_subjects = LENGTH(Svoxel);

  PROTECT(Soutput=allocVector(REALSXP, *n_groups));
  sum_voxel = REAL(Soutput);


  /* init variables */
  for(i=0; i < *n_groups; i++) {
    sum_voxel[i] = 0;
  }

  for(i=0; i < n_subjects; i++) {
    sum_voxel[(int) groupings[i]] += voxel[i];
  }

  UNPROTECT(1);
  return(Soutput);
}

SEXP voxel_mean(SEXP Svoxel, SEXP Sn_groups, SEXP Sgroupings) {
  double *voxel, *groupings, *n_groups, *means;
  double *sum_voxel;
  int *subjects_per_group;
  int n_subjects, i;
  SEXP Soutput;


  voxel = REAL(Svoxel);
  n_groups = REAL(Sn_groups);
  groupings = REAL(Sgroupings);

  n_subjects = LENGTH(Svoxel);

  PROTECT(Soutput=allocVector(REALSXP, *n_groups));
  means = REAL(Soutput);

  subjects_per_group = malloc(sizeof(int) * *n_groups);
  sum_voxel = malloc(sizeof(double) * *n_groups);

  /* init variables */
  for(i=0; i < *n_groups; i++) {
    subjects_per_group[i] = 0;
    sum_voxel[i] = 0;
  }

  for(i=0; i < n_subjects; i++) {
    subjects_per_group[(int) groupings[i]]++;
    sum_voxel[(int) groupings[i]] += voxel[i];
  }

  for(i=0; i< *n_groups; i++) {
    means[i] = sum_voxel[i] / subjects_per_group[i];
  }

  free(sum_voxel);
  free(subjects_per_group);
  UNPROTECT(1);
  return(Soutput);
}

SEXP voxel_var(SEXP Svoxel, SEXP Sn_groups, SEXP Sgroupings) {
  double *voxel, *groupings, *n_groups, *means;
  double *var, *sum_voxel;
  int *subjects_per_group;
  int n_subjects, i;
  SEXP Soutput;
  SEXP Smean;


  voxel = REAL(Svoxel);
  n_groups = REAL(Sn_groups);
  groupings = REAL(Sgroupings);

  n_subjects = LENGTH(Svoxel);

  PROTECT(Soutput=allocVector(REALSXP, *n_groups));
  var = REAL(Soutput);

  subjects_per_group = malloc(sizeof(int) * *n_groups);
  sum_voxel = malloc(sizeof(double) * *n_groups);

  Smean = voxel_mean(Svoxel, Sn_groups, Sgroupings);
  means = REAL(Smean);

  /* init variables */
  for(i=0; i < *n_groups; i++) {
    subjects_per_group[i] = 0;
    sum_voxel[i] = 0;
  }

  for(i=0; i < n_subjects; i++) {
    subjects_per_group[(int) groupings[i]]++;
    sum_voxel[(int) groupings[i]] 
      += pow(voxel[i] - means[(int) groupings[i]], 2);
  }



  for(i=0; i< *n_groups; i++) {
    var[i] = sum_voxel[i] / (subjects_per_group[i]-1);
  }
  free(sum_voxel);
  free(subjects_per_group);
  UNPROTECT(1);
  return(Soutput);
}


SEXP voxel_correlation(SEXP Sx, SEXP Sy) {
  double *x, *y,*r;
  double sum_x, sum_y, sum_xy, sum_x2, sum_y2, numerator, denominator, n;
  int i;
  SEXP output;

  n = LENGTH(Sy);
  PROTECT(output=allocVector(REALSXP, 1));
  
  r = REAL(output);
  x = REAL(Sx);
  y = REAL(Sy);

  /* compute sums for x and y */
  sum_x = 0;
  sum_y = 0;
  sum_xy = 0;
  sum_y2 = 0;
  sum_x2 = 0;
  for (i=0; i < n; i++) {
    sum_x += x[i];
    sum_y += y[i];
    sum_xy += x[i] * y[i];
    sum_x2 += pow(x[i], 2);
    sum_y2 += pow(y[i], 2);
  }

  numerator = n * sum_xy - sum_x * sum_y;
  denominator = sqrt( (n * sum_x2 - pow(sum_x, 2)) * 
		      (n * sum_y2 - pow(sum_y, 2)) );
  *r = numerator / denominator;
  UNPROTECT(1);
  return(output);
}

SEXP voxel_lm_file(SEXP Sy, SEXP Sx,int n,int p,double *coefficients, 
        double *residuals, double *effects, 
        double *work, 
        double *qraux, double *v, int *pivot, double *se, double *t) {
 
  double tol, rss, resvar, mss, mean_fitted, sum_fitted;
  double *x, *y, *xoutput;
  int ny,ny1, rank, i, j, rdf, index;
  SEXP new_x, output;
  int nprot = 0;


  ny = ncols(Sy);
  ny1 = nrows(Sy);


  //Rprintf("n %d p %d\n", n,p);

  // the output will contain:
  // 
  // f-statistic
  // p * t-statistic
  // r-squared
  //
  // which is a total of p+2 values
  PROTECT(output=allocVector(REALSXP, p+2)); nprot++;
  xoutput = REAL(output);

  /* since x (the model matrix, input variable Sx) is destroyed in dqrls, 
     create a copy here */
  PROTECT(new_x=allocMatrix(REALSXP, n, p)); nprot++;
  for (i=0; i < n; i++) {
    for (j=0; j < p; j++) {
      REAL(new_x)[i+j*n] = REAL(Sx)[i+j*n];
      //Rprintf("new_x %f\n", REAL(new_x)[i+j*n]);
    }
  }

  x = REAL(new_x);
  y = REAL(Sy);

  //Rprintf("coly %d rowy %d\n", ny,ny1);
  rank = 1;
  tol = 1e-07;

  // compute the least squares solution:
  F77_CALL(dqrls)(x, &n, &p, y, &ny, &tol, coefficients, residuals, effects,
        &rank, pivot, qraux, work);

  // Calculate the f-statistic first
  rss = 0; // residual sum of squares
  rdf = 0; // residual degrees of freedom
  mss = 0; // (fitted - mean fitted sum) of squares
  sum_fitted = 0;
  for (i=0; i < n; i++) {
    rss += pow(residuals[i], 2);
    sum_fitted += (y[i] - residuals[i]);
  }

  mean_fitted = sum_fitted / n;
  for (i=0; i < n; i++) {
    mss += pow((y[i] - residuals[i]) - mean_fitted, 2);
  }
  
  rdf = n - p;

  //Rprintf("rss %f p %d resvar %f %\n", rss,p,resvar);
  resvar = rss/rdf;
  //Rprintf("mss %f p %d resvar %f %\n", mss,p,resvar);

  /* first output is the f-stat of the whole model */

  
  xoutput[0] = (mss/(p - 1))/resvar;
  
  // DPOTRI - compute the inverse of a real symmetric positive
  // definite matrix A using the Cholesky factorization A =
  // U**T*U or A = L*L**T computed by DPOTRF
  int info;
  F77_CALL(dpotri)("Upper", &p, x, &n, &info);

  if (info != 0) {
    UNPROTECT(nprot);
    if (info > 0) {
      error(("element (%d, %d) is zero, so the inverse cannot be computed"),info, info);
    }
    error(("argument %d of Lapack routine %s had invalid value"), -info, "dpotri");
  }

  // on to the t-statistics for the intercept and other terms
  for (i=0; i < p; i++) {
    index = (i * n) + i;
    se[i] = sqrt(x[index] *resvar);
    xoutput[i+1] = coefficients[i] / se[i];
  }

  // last, but not least, the r-squared:
  xoutput[p+1] = mss / (mss + rss);
  
  UNPROTECT(nprot);
  return(output);
}


SEXP voxel_lm(SEXP Sy, SEXP Sx, double *coefficients,
        double *residuals, double *effects,
        double *work,
        double *qraux, double *v, int *pivot, double *se, double *t) {

  double tol, rss, resvar, mss, mean_fitted, sum_fitted;
  double *x, *y, *xoutput;
  int n, p, ny, rank, i, j, rdf, index;
  SEXP new_x, output;
  int nprot = 0;

  n = nrows(Sx);
  p = ncols(Sx);

  // the output will contain:
  //
  // f-statistic
  // p * t-statistic
  // r-squared
  //
  // which is a total of p+2 values
  PROTECT(output=allocVector(REALSXP, p+2)); nprot++;
  xoutput = REAL(output);

  /* since x (the model matrix, input variable Sx) is destroyed in dqrls,
create a copy here */
  PROTECT(new_x=allocMatrix(REALSXP, n, p)); nprot++;
  for (i=0; i < n; i++) {
    for (j=0; j < p; j++) {
      REAL(new_x)[i+j*n] = REAL(Sx)[i+j*n];
    }
  }

  x = REAL(new_x);
  y = REAL(Sy);
  ny = ncols(Sy);
  rank = 1;
  tol = 1e-07;

  // compute the least squares solution:
  F77_CALL(dqrls)(x, &n, &p, y, &ny, &tol, coefficients, residuals, effects,
        &rank, pivot, qraux, work);

  // Calculate the f-statistic first
  rss = 0; // residual sum of squares
  rdf = 0; // residual degrees of freedom
  mss = 0; // (fitted - mean fitted sum) of squares
  sum_fitted = 0;
  for (i=0; i < n; i++) {
    rss += pow(residuals[i], 2);
    sum_fitted += (y[i] - residuals[i]);
  }
  mean_fitted = sum_fitted / n;
  for (i=0; i < n; i++) {
    mss += pow((y[i] - residuals[i]) - mean_fitted, 2);
  }
  rdf = n - p;
  resvar = rss/rdf;
  /* first output is the f-stat of the whole model */
  xoutput[0] = (mss/(p - 1))/resvar;

  // DPOTRI - compute the inverse of a real symmetric positive
  // definite matrix A using the Cholesky factorization A =
  // U**T*U or A = L*L**T computed by DPOTRF
  int info;
  F77_CALL(dpotri)("Upper", &p, x, &n, &info);

  if (info != 0) {
    UNPROTECT(nprot);
    if (info > 0) {
      error(("element (%d, %d) is zero, so the inverse cannot be computed"),info, info);
    }
    error(("argument %d of Lapack routine %s had invalid value"), -info, "dpotri");
  }

  // on to the t-statistics for the intercept and other terms
  for (i=0; i < p; i++) {
    index = (i * n) + i;
    se[i] = sqrt(x[index] *resvar);
    xoutput[i+1] = coefficients[i] / se[i];
  }

  // last, but not least, the r-squared:
  xoutput[p+1] = mss / (mss + rss);
  
  UNPROTECT(nprot);
  return(output);
}


/*
SEXP voxel_anova(SEXP Sy, SEXP Sx, SEXP asgn,
		 double *coefficients, 
		 double *residuals, 
		 double *effects, 
		 double *work, 
		 double *qraux, 
		 double *v, 
		 int *pivot, 
		 double *se, 
		 double *t,
		 double *comp,
		 double *ss,
		 int *df) {

  int p, n, maxasgn, dfr, i;
  int *xasgn;
  double ssr;
  double *xf_sexp;
  SEXP f_sexp;

  n = nrows(Sx);
  p = ncols(Sx);


  //Rprintf("N: %d P: %d\n", n,p);
  xasgn = INTEGER(asgn);

  maxasgn = 0;
  for (i=0; i < p; i++) {
    if (xasgn[i] > maxasgn) {
      maxasgn = (int) xasgn[i];
    }
  }
  maxasgn++;

  PROTECT(f_sexp = allocVector(REALSXP, maxasgn-1));
  xf_sexp = REAL(f_sexp);
  
  voxel_lm(Sy, Sx, coefficients, residuals, effects,
	   work, qraux, v, pivot, se, t);

  // dfr: degrees of freedom of the residuals
  dfr = n - p;
  for (i=0; i < p; i++) {
    comp[i] = pow(effects[i], 2);
  }

  // set sum of squares and degrees of freedom to zero
  for (i=0; i < maxasgn; i++) {
    ss[i] = 0.0;
    df[i] = 0;
  }

  // compute sums of squares
  for (i=0; i < p; i++) {
    ss[(int)xasgn[i]] += comp[i];
    df[xasgn[i]]++;
    //Rprintf("%d %d %f %f\n", i, xasgn[i], comp[i], ss[xasgn[i]]);
  }
  
  // ssr: sums of squares of residuals
  ssr = 0.0;
  for (i=0; i < n; i++) {
    ssr += pow(residuals[i], 2);
  }
  // compute the f-statistics
  for (i=1; i<maxasgn; i++) {
    xf_sexp[i-1] = (ss[i] / df[i]) / (ssr/dfr);
    //Rprintf("F value: %f\n", xf_sexp[i-1]);
  }

  //Rprintf("effects: %f\n", effects[0]);
  UNPROTECT(1);
  return(f_sexp);
}

SEXP test_voxel_anova(SEXP Sy, SEXP Sx, SEXP asgn) {
  int                result;
  mihandle_t         *hvol, hmask;
  char               *method_name;
  int                i, v0, v1, v2, output_index, buffer_index;
  unsigned long      start[3], count[3];
  unsigned long      location[3];
  int                num_files;
  double             *xn_groups;
  double             *xbuffer, *xoutput, **full_buffer, *xhave_mask, *xn;
  double             *mask_buffer;
  double             *groupings;
  midimhandle_t      dimensions[3];
  unsigned int       sizes[3];
  SEXP               output, buffer, R_fcall, t_sexp, n_groups, f_sexp;

  double             *y, *x, *coefficients, *residuals, *effects;
  double             *diag, *se, *t, *work, *qraux, *v, *comp, *xf_sexp, *ss;
  double             tol, rss, resvar, ssr;
  int                n, p, ny, rank, j, rdf, index, dfr, maxasgn;
  int                *pivot, *xasgn, *df;


  n = nrows(Sx);
  p = ncols(Sx);

  coefficients = malloc(sizeof(double) * p);
  residuals = malloc(sizeof(double) * n);
  effects = malloc(sizeof(double) * n);
  pivot = malloc(sizeof(int) * p);
  work = malloc(sizeof(double) * (2*p));
  qraux = malloc(sizeof(double) * p);
  v = malloc(sizeof(double) * p * p);
  diag = malloc(sizeof(double) * p);
  se = malloc(sizeof(double) * p);
  t = malloc(sizeof(double) * p);
  comp = malloc(sizeof(double) * p);
  ss = malloc(sizeof(double) * maxasgn);
  df = malloc(sizeof(int) * maxasgn);
  xasgn = INTEGER(asgn);

  maxasgn = 0;
  for (i=0; i < p; i++) {
    if (xasgn[i] > maxasgn) {
      maxasgn = (int) xasgn[i];
    }
  }
  maxasgn++;

  //PROTECT(f_sexp = allocVector(REALSXP, maxasgn-1));
  f_sexp = voxel_anova(Sy, Sx, asgn,
		       coefficients, residuals,
		       effects, work, qraux, v, pivot,
		       se, t, comp, ss, df);
  //UNPROTECT(1);

  free(coefficients);
  free(residuals);
  free(effects);
  free(pivot);
  free(work);
  free(qraux);
  free(v);
  free(diag);
  free(se);
  free(t);
  free(comp);
  free(ss);
  free(df);

  return(f_sexp);
}
*/
  
  
/* minc2_model: run one of a set of modelling function at every voxel
 * filenames: character list of minc2 volumes.
 * Sx: variable that each particular function will work on. The model matrix 
 *     for method "lm", for example, or the function string for method "eval"
 * asgn: assignments for factors - only used by anova.
 * have_mask: a double of either 0 or 1 depending on whether a mask should
 *            be used.
 * mask: a string containing the mask filename.
 * mask_value: the value inside the mask at which the function is to be evaled
 * rho: the environment - only used by method "eval"
 * nresults: the number of columns in the result - only used by method "eval"
 * method: a string containing one of "t-test", "wilcoxon", or "correlation"
 */
SEXP minc2_model(SEXP filenames, SEXP Sx, SEXP asgn,
		 SEXP have_mask, SEXP mask, SEXP mask_lower_value,
		 SEXP mask_upper_value, SEXP rho, SEXP nresults, SEXP method) {
  int                result;
  mihandle_t         *hvol, hmask;
  char               *method_name;
  int                i, v0, v1, v2, output_index, buffer_index;
  unsigned long      start[3], count[3];
  unsigned long      location[3];
  int                num_files;
  double             *xn_groups;
  double             *xbuffer, *xoutput, **full_buffer, *xhave_mask;
  double             *xmask_lower_value;
  double             *xmask_upper_value;
  double             *mask_buffer;
  double             *groupings;
  midimhandle_t      dimensions[3];
  unsigned int       sizes[3];
  SEXP               output, buffer, t_sexp, n_groups;
  /* stuff for linear models only */
  double             *coefficients, *residuals, *effects; 
  double             *diag, *se, *t, *work, *qraux, *v, *ss, *comp;
  int                n, p, maxasgn;
  int                *pivot, *xasgn, *df;

  num_files = LENGTH(filenames);

  /* get the method that should be used at each voxel */
  method_name = CHAR(STRING_ELT(method, 0));
  Rprintf("Method: %s\n", method_name);

  /* allocate memory for the volume handles */
  hvol = malloc(num_files * sizeof(mihandle_t));

  Rprintf("Number of volumes: %i\n", num_files);

  /* open the mask - if so desired */
  xhave_mask = REAL(have_mask);
  if (xhave_mask[0] == 1) {
    result = miopen_volume(CHAR(STRING_ELT(mask, 0)),
			   MI2_OPEN_READ, &hmask);
    if (result != MI_NOERROR) {
      error("Error opening mask: %s.\n", CHAR(STRING_ELT(mask, 0)));
    }
  }
  
  /* get the value at which the mask is to be evaluated */
  xmask_lower_value = REAL(mask_lower_value);
  xmask_upper_value = REAL(mask_upper_value);

  /* open each volume */
  for(i=0; i < num_files; i++) {
    result = miopen_volume(CHAR(STRING_ELT(filenames, i)),
      MI2_OPEN_READ, &hvol[i]);
    if (result != MI_NOERROR) {
      error("Error opening input file: %s.\n", CHAR(STRING_ELT(filenames,i)));
    }
  }

  /* get the file dimensions and their sizes - assume they are the same*/
  miget_volume_dimensions( hvol[0], MI_DIMCLASS_SPATIAL,
			   MI_DIMATTR_ALL, MI_DIMORDER_FILE,
			   3, dimensions);
  result = miget_dimension_sizes( dimensions, 3, sizes );
  Rprintf("Volume sizes: %i %i %i\n", sizes[0], sizes[1], sizes[2]);


  /* allocate the local buffer that will be passed to the function */
  PROTECT(buffer=allocVector(REALSXP, num_files));
  xbuffer = REAL(buffer); 

  /* allocate stuff for means and standard deviations */
  if (strcmp(method_name, "mean") == 0 ||
      strcmp(method_name, "sum") == 0 ||
      strcmp(method_name, "var") == 0) {
    PROTECT(n_groups=allocVector(REALSXP, 1));
    xn_groups = REAL(n_groups);
    groupings = REAL(Sx);

    xn_groups[0] = 0;

    /* determine the number of groups. Here we assume that the input
       will contain integers corresponding to the group number, so we need
       but take the max of that to get the number of groups. Note that we
       have to add 1 to the number of groups, as the first group will have 
       a value of 0 */
    for(i=0; i < LENGTH(Sx); i++) {
      if (groupings[i] > xn_groups[0]) {
	xn_groups[0]++;
      }
    }
    xn_groups[0]++;

    Rprintf("N GROUPS: %f\n", xn_groups[0]);
    PROTECT(t_sexp = allocVector(REALSXP, xn_groups[0]));
    PROTECT(output=allocMatrix(REALSXP, (sizes[0] * sizes[1] * sizes[2]), 
			       xn_groups[0]));

  }
  /* allocate stuff for evaluation arbitrary functions */
  else if (strcmp(method_name, "eval") == 0) {
    xn_groups = REAL(nresults);
    PROTECT(t_sexp = allocVector(REALSXP, xn_groups[0]));
    PROTECT(output=allocMatrix(REALSXP, (sizes[0] * sizes[1] * sizes[2]),
			       xn_groups[0]));
  }
  /* allocate stuff necessary for fitting linear models */
  else if (strcmp(method_name, "lm") == 0) {
    n = nrows(Sx);
    p = ncols(Sx);
    
    coefficients = malloc(sizeof(double) * p);
    residuals = malloc(sizeof(double) * n);
    effects = malloc(sizeof(double) * n);
    pivot = malloc(sizeof(int) * p);
    work = malloc(sizeof(double) * (2*p));
    qraux = malloc(sizeof(double) * p);
    v = malloc(sizeof(double) * p * p);
    diag = malloc(sizeof(double) * p);
    se = malloc(sizeof(double) * p);
    t = malloc(sizeof(double) * p);
    
    Rprintf("N: %d P: %d\n", n,p);

    PROTECT(t_sexp = allocVector(REALSXP, p + 2));

    /* allocate the output buffer */
    PROTECT(output=allocMatrix(REALSXP, (sizes[0] * sizes[1] * sizes[2]), 2*p + 2));

  }
  else if (strcmp(method_name, "anova") == 0) {
    n = nrows(Sx);
    p = ncols(Sx);
    xasgn = INTEGER(asgn);

    maxasgn = 0;
    for (i=0; i < p; i++) {
      if (xasgn[i] > maxasgn) {
	maxasgn = (int) xasgn[i];
      }
    }
    maxasgn++;

    coefficients = malloc(sizeof(double) * p);
    residuals = malloc(sizeof(double) * n);
    effects = malloc(sizeof(double) * n);
    pivot = malloc(sizeof(int) * p);
    work = malloc(sizeof(double) * (2*p));
    qraux = malloc(sizeof(double) * p);
    v = malloc(sizeof(double) * p * p);
    diag = malloc(sizeof(double) * p);
    se = malloc(sizeof(double) * p);
    t = malloc(sizeof(double) * p);

    comp = malloc(sizeof(double) * p);
    ss = malloc(sizeof(double) * maxasgn);
    df = malloc(sizeof(int) * maxasgn);
    
    Rprintf("N: %d P: %d\n", n,p);

    PROTECT(t_sexp = allocVector(REALSXP, maxasgn-1));

    /* allocate the output buffer */
    PROTECT(output=allocMatrix(REALSXP, (sizes[0] * sizes[1] * sizes[2]), 
			       maxasgn-1));

  }
  else {    
    /* allocate the output buffer */
    PROTECT(output=allocVector(REALSXP, (sizes[0] * sizes[1] * sizes[2])));
  }
  xoutput = REAL(output);
    
  //PROTECT(R_fcall = lang2(fn, R_NilValue));


  /* allocate first dimension of the buffer */
  full_buffer = malloc(num_files * sizeof(double));

  /* allocate second dimension of the buffer 
     - big enough to hold one slice per subject at a time */
  for (i=0; i < num_files; i++) {
    full_buffer[i] = malloc(sizes[1] * sizes[2] * sizeof(double));
  }
  
  /* allocate buffer for mask - if necessary */
  if (xhave_mask[0] == 1) {
    mask_buffer = malloc(sizes[1] * sizes[2] * sizeof(double));
  }
	
  /* set start and count. start[0] will change during the loop */
  start[0] = 0; start[1] = 0; start[2] = 0;
  count[0] = 1; count[1] = sizes[1]; count[2] = sizes[2];

  /* loop across all files and voxels */
  Rprintf("In slice \n");
  for (v0=0; v0 < sizes[0]; v0++) {
    start[0] = v0;
    for (i=0; i < num_files; i++) {
      if (miget_real_value_hyperslab(hvol[i], 
				     MI_TYPE_DOUBLE, 
				     (unsigned long *) start, 
				     (unsigned long *) count, 
				     full_buffer[i]) )
	error("Error opening buffer.\n");
    }
    /* get mask - if desired */
    if (xhave_mask[0] == 1) {
      if (miget_real_value_hyperslab(hmask, 
				     MI_TYPE_DOUBLE, 
				     (unsigned long *) start, 
				     (unsigned long *) count, 
				     mask_buffer) )
	error("Error opening mask buffer.\n");
    }

    Rprintf(" %d ", v0);
    for (v1=0; v1 < sizes[1]; v1++) {
      for (v2=0; v2 < sizes[2]; v2++) {
	output_index = v0*sizes[1]*sizes[2]+v1*sizes[2]+v2;
	buffer_index = sizes[2] * v1 + v2;

	/* only perform operation if not masked */
	if(xhave_mask[0] == 0 
	   || (xhave_mask[0] == 1 && 
	       mask_buffer[buffer_index] > xmask_lower_value[0] -0.5 &&
	       mask_buffer[buffer_index] < xmask_upper_value[0] +0.5)) {
	
	  for (i=0; i < num_files; i++) {
	    location[0] = v0;
	    location[1] = v1;
	    location[2] = v2;

	    xbuffer[i] = full_buffer[i][buffer_index];
	    
	    //Rprintf("V%i: %f\n", i, full_buffer[i][index]);

	  }

	  /* compute either a t test of wilcoxon rank sum test */
	  if (strcmp(method_name, "t-test") == 0) {
	    xoutput[output_index] = REAL(t_test(buffer, Sx))[0]; 
	  }
	  else if (strcmp(method_name, "paired-t-test") == 0) {
	    xoutput[output_index] = REAL(paired_t_test(buffer, Sx))[0];
	  }
	  else if (strcmp(method_name, "wilcoxon") == 0) {
	    xoutput[output_index] = 
	      REAL(wilcoxon_rank_test(buffer, Sx))[0];
	  }
	  else if (strcmp(method_name, "correlation") == 0) {
	    xoutput[output_index] = 
	      REAL(voxel_correlation(buffer, Sx))[0];
	  }
	  else if (strcmp(method_name, "mean") == 0) {
	    t_sexp = voxel_mean(buffer, n_groups, Sx);
	    for(i=0; i < xn_groups[0]; i++) {
	      xoutput[output_index + i * (sizes[0]*sizes[1]*sizes[2])] 
		= REAL(t_sexp)[i];
	    }
	  }
	  else if (strcmp(method_name, "sum") == 0) {
	    t_sexp = voxel_sum(buffer, n_groups, Sx);
	    for(i=0; i < xn_groups[0]; i++) {
	      xoutput[output_index + i * (sizes[0]*sizes[1]*sizes[2])] 
		= REAL(t_sexp)[i];
	    }
	  }
	  else if (strcmp(method_name, "var") == 0) {
	    t_sexp = voxel_var(buffer, n_groups, Sx);
	    for(i=0; i < xn_groups[0]; i++) {
	      xoutput[output_index + i * (sizes[0]*sizes[1]*sizes[2])] 
		= REAL(t_sexp)[i];
	    }
	  }
	  else if (strcmp(method_name, "eval") == 0) {
	    /* install the variable "x" into environment */
	    defineVar(install("x"), buffer, rho);
	    t_sexp = eval(Sx, rho);
	    for(i=0; i < xn_groups[0]; i++) {
	      xoutput[output_index + i * (sizes[0]*sizes[1]*sizes[2])] 
		= REAL(t_sexp)[i];
	    }
	  }
    else if (strcmp(method_name, "lm") == 0) {
    t_sexp = voxel_lm(buffer, Sx, coefficients, residuals, effects,
            work, qraux, v, pivot, se, t);
      
      // most sensible output format (?): fist the full model measurements,
      // then the individual measurement in the same order as summary.lm
      // gives them:
      //
      // f-statistic
      // r-squared
      // betas
      // t-stats
      //
      
      // f-statistic
      xoutput[output_index] = REAL(t_sexp)[0];
      
      // r-squared (last value from voxel_lm call: p+2 (stating at 0, so p+1))
      xoutput[output_index + 1 * (sizes[0]*sizes[1]*sizes[2])] = REAL(t_sexp)[p + 1];
      
      // the betas/coefficients:
      for (int k = 2; k < (p + 2); k++) {
        xoutput[output_index + k * (sizes[0]*sizes[1]*sizes[2])] = coefficients[k - 2];
      }
      
      // t-stats
      for(int k = 1; k < p + 1; k++) {
        xoutput[output_index + (k + p + 1) * (sizes[0]*sizes[1]*sizes[2])] = REAL(t_sexp)[k];
      }
  }
	  /*
	  else if (strcmp(method_name, "anova") == 0) {
	    t_sexp = voxel_anova(buffer, Sx, asgn,
				 coefficients, residuals,
				 effects, work, qraux, v, pivot,
				 se, t, comp, ss, df);
	    for(i=0; i < maxasgn-1; i++) {
	      xoutput[output_index + i * (sizes[0]*sizes[1]*sizes[2])] 
		      = REAL(t_sexp)[i];
	    }
	  }
	  */
	}
	else {
	  xoutput[output_index] = 0;
	}
      }
    }
  }
  Rprintf("\nDone\n");

  /* free memory */
  for (i=0; i<num_files; i++) {
    miclose_volume(hvol[i]);
    free(full_buffer[i]);
  }
  if (strcmp(method_name, "lm")==0) {
    free(full_buffer);
    free(coefficients);
    free(residuals);
    free(effects);
    free(pivot);
    free(work);
    free(qraux);
    free(v);
    free(diag);
    free(se);
    free(t);
    UNPROTECT(3);
  }
  else if (strcmp(method_name, "anova")==0) {
    free(full_buffer);
    free(coefficients);
    free(residuals);
    free(effects);
    free(pivot);
    free(work);
    free(qraux);
    free(v);
    free(diag);
    free(se);
    free(t);
    free(comp);
    free(ss);
    free(df);
    UNPROTECT(3);
  }

  else if (strcmp(method_name, "mean") == 0 ||
	   strcmp(method_name, "var") == 0 ||
	   strcmp(method_name, "sum") == 0) {
    UNPROTECT(4);
  }
  else if (strcmp(method_name, "eval") == 0) {
    UNPROTECT(3);
  }
  else {
    UNPROTECT(2);
  }

  /* return the results */
  return(output);
}

/* minc2_model: run one of a set of modelling function at every voxel
 * filenames: character list of minc2 volumes.
 * Sx: variable that each particular function will work on. The model matrix 
 *     for method "lm", for example, or the function string for method "eval"
 * asgn: assignments for factors - only used by anova.
 * have_mask: a double of either 0 or 1 depending on whether a mask should
 *            be used.
 * mask: a string containing the mask filename.
 * mask_value: the value inside the mask at which the function is to be evaled
 * rho: the environment - only used by method "eval"
 * nresults: the number of columns in the result - only used by method "eval"
 * method: a string containing one of "t-test", "wilcoxon", or "correlation"
 */
SEXP minc2_model_file(SEXP filenames_left,SEXP filenames_right, SEXP mmatrix, SEXP asgn,
		 SEXP have_mask, SEXP mask, SEXP mask_lower_value,
		 SEXP mask_upper_value, SEXP rho, SEXP nresults, SEXP method) {
  int                result;
  mihandle_t         *hvol, *hvol_left,*hvol_right,hmask;
  char               *method_name;
  int                i, v0, v1, v2, output_index, buffer_index;
  unsigned long      start[3], count[3];
  unsigned long      location[3];
  int                num_files_left,num_files_right;
  double             *xn_groups;
  double             *xbuffer,*ybuffer, *xoutput, **full_buffer_left,**full_buffer_right, *xhave_mask;
  double             *xmask_lower_value;
  double             *xmask_upper_value;
  double             *mask_buffer;
  double             *groupings;
  midimhandle_t      dimensions[3];
  unsigned int       sizes[3];
  SEXP               output, buffer,buffer1, t_sexp, n_groups;
  /* stuff for linear models only */
  double             *coefficients, *residuals, *effects; 
  double             *diag, *se, *t, *work, *qraux, *v, *ss, *comp;
  int                n, p, maxasgn,mmatrix_rows,mmatrix_cols;
  int                *pivot, *xasgn, *df;
  double 	     *pMmatrix;
  num_files_left = LENGTH(filenames_left);

  /* allocate memory for the volume handles */
  hvol_left = malloc(num_files_left * sizeof(mihandle_t));
  hvol_right = malloc(num_files_left * sizeof(mihandle_t));
  Rprintf("Number of volumes: %i\n", num_files_left);

  /* open the mask - if so desired */
  xhave_mask = REAL(have_mask);
  if (xhave_mask[0] == 1) {
    result = miopen_volume(CHAR(STRING_ELT(mask, 0)),
			   MI2_OPEN_READ, &hmask);
    if (result != MI_NOERROR) {
      error("Error opening mask: %s.\n", CHAR(STRING_ELT(mask, 0)));
    }
  }
  
  /* get the value at which the mask is to be evaluated */
  xmask_lower_value = REAL(mask_lower_value);
  xmask_upper_value = REAL(mask_upper_value);

  /* open each volume */
  for(i=0; i < num_files_left; i++) {
    result = miopen_volume(CHAR(STRING_ELT(filenames_left, i)),
      MI2_OPEN_READ, &hvol_left[i]);
    if (result != MI_NOERROR) {
      error("Error opening input file: %s.\n", CHAR(STRING_ELT(filenames_left,i)));
    }
  }
  for(i=0; i < num_files_left; i++) {
    result = miopen_volume(CHAR(STRING_ELT(filenames_right, i)),
      MI2_OPEN_READ, &hvol_right[i]);
    if (result != MI_NOERROR) {
      error("Error opening input file: %s.\n", CHAR(STRING_ELT(filenames_right,i)));
    }
  }

  /* get the file dimensions and their sizes - assume they are the same*/
  miget_volume_dimensions( hvol_left[0], MI_DIMCLASS_SPATIAL,
			   MI_DIMATTR_ALL, MI_DIMORDER_FILE,
			   3, dimensions);
  result = miget_dimension_sizes( dimensions, 3, sizes );
  Rprintf("Volume sizes: %i %i %i\n", sizes[0], sizes[1], sizes[2]);


  // determine numbers of rows, columns, and vertices

  // Case 1: There is no static part
  if(isLogical(mmatrix)) 
	// For now, maximum allowed dynamic parts is 1 so set p = 2 (1 for intercept)
	p = 2;
  // Case 2: There is a static part
  else {
        pMmatrix = REAL(mmatrix);
  	mmatrix_cols = ncols(mmatrix);
  	mmatrix_rows = nrows(mmatrix);
        p = mmatrix_cols + 1;
        Rprintf("mmatrix cols: %d mmatrix rows: %d\n", mmatrix_cols,mmatrix_rows );
  }
  n = num_files_left;
  /* allocate the local buffer that will be passed to the function */
  PROTECT(buffer=allocVector(REALSXP, num_files_left));
  PROTECT(buffer1=allocVector(REALSXP, num_files_left*p));
  xbuffer = REAL(buffer); 
  ybuffer=REAL(buffer1);


   
  coefficients = malloc(sizeof(double) * p);
  residuals = malloc(sizeof(double) * n);
  effects = malloc(sizeof(double) * n);
  pivot = malloc(sizeof(int) * p);
  work = malloc(sizeof(double) * (2*p));
  qraux = malloc(sizeof(double) * p);
  v = malloc(sizeof(double) * p * p);
  diag = malloc(sizeof(double) * p);
  se = malloc(sizeof(double) * p);
  t = malloc(sizeof(double) * p);
    
  Rprintf("N: %d P: %d\n", n,p);

  PROTECT(t_sexp = allocVector(REALSXP, p + 2));
 


  /* allocate the output buffer */
  PROTECT(output=allocMatrix(REALSXP, (sizes[0] * sizes[1] * sizes[2]), 2*p + 2));

  
  xoutput = REAL(output);
    
  //PROTECT(R_fcall = lang2(fn, R_NilValue));


  /* allocate first dimension of the buffer */
  full_buffer_left = malloc(num_files_left * sizeof(double));
  full_buffer_right = malloc(num_files_left * sizeof(double));
  /* allocate second dimension of the buffer 
     - big enough to hold one slice per subject at a time */
  for (i=0; i < num_files_left; i++) {
    full_buffer_left[i] = malloc(sizes[1] * sizes[2] * sizeof(double));
    full_buffer_right[i] = malloc(sizes[1] * sizes[2] * sizeof(double));
  }
  
  /* allocate buffer for mask - if necessary */
  if (xhave_mask[0] == 1) {
    mask_buffer = malloc(sizes[1] * sizes[2] * sizeof(double));
  }
	
  /* set start and count. start[0] will change during the loop */
  start[0] = 0; start[1] = 0; start[2] = 0;
  count[0] = 1; count[1] = sizes[1]; count[2] = sizes[2];

  /* loop across all files and voxels */


  Rprintf("In slice \n");
  for (v0=0; v0 < sizes[0]; v0++) {
    start[0] = v0;
    for (i=0; i < num_files_left; i++) {
      if (miget_real_value_hyperslab(hvol_left[i], 
				     MI_TYPE_DOUBLE, 
				     (unsigned long *) start, 
				     (unsigned long *) count, 
				     full_buffer_left[i]) )
	error("Error opening buffer.\n");
      if (miget_real_value_hyperslab(hvol_right[i], 
				     MI_TYPE_DOUBLE, 
				     (unsigned long *) start, 
				     (unsigned long *) count, 
				     full_buffer_right[i]) )
	error("Error opening buffer.\n");
    }
   
    if (xhave_mask[0] == 1) {
      if (miget_real_value_hyperslab(hmask, 
				     MI_TYPE_DOUBLE, 
				     (unsigned long *) start, 
				     (unsigned long *) count, 
				     mask_buffer) )
	error("Error opening mask buffer.\n");
    }

    Rprintf(" %d ", v0);
    for (v1=0; v1 < sizes[1]; v1++) {
      for (v2=0; v2 < sizes[2]; v2++) {
	output_index = v0*sizes[1]*sizes[2]+v1*sizes[2]+v2;
	buffer_index = sizes[2] * v1 + v2;

	if(xhave_mask[0] == 0 
	   || (xhave_mask[0] == 1 && 
	       mask_buffer[buffer_index] > xmask_lower_value[0] -0.5 &&
	       mask_buffer[buffer_index] < xmask_upper_value[0] +0.5)) {
	
	  for (i=0; i < num_files_left; i++) {
	    location[0] = v0;
	    location[1] = v1;
	    location[2] = v2;

	    xbuffer[i] = full_buffer_left[i][buffer_index];
	    }

	    if(isLogical(mmatrix))  { 
		    // fill y buffer
		    // Intercept
		    for (int j=0; j < n; j++) {
			ybuffer[j] = 1.0;
			//Rprintf("ybuffer %f index %d\n", ybuffer[j] ,j);
		    }    
		    // Current Vertex Data
		    for (int j=0; j<n; j++) {
		      ybuffer[j+n] = full_buffer_right[j][buffer_index];
		      //Rprintf("ybuffer %f index %d\n", ybuffer[j+n],j+n);
		    }    
	    }
	    else {
		// Fill with static part
	    	for (int j=0; j < mmatrix_cols*mmatrix_rows; j++) {
			ybuffer[j] = pMmatrix[j];
		        //Rprintf("mmatrix %f index %d\n", pMmatrix[j],j);
		    }  
		// Fill with dynamic part  
	    	for (int j=0; j < n ; j++) {
			ybuffer[j+mmatrix_cols*mmatrix_rows] = full_buffer_right[j][buffer_index];
		        //Rprintf("mmatrix %f index %d\n", ydata[i+nVertices*j],j+mmatrix_cols*mmatrix_rows);
		    }   
	    }


	    //Rprintf("V%i: %f\n", i, full_buffer[i][index]);


    t_sexp = voxel_lm_file(buffer, buffer1,n,p ,coefficients, residuals, effects,
          work, qraux, v, pivot, se, t);
      
      // most sensible output format (?): fist the full model measurements,
      // then the individual measurement in the same order as summary.lm
      // gives them:
      //
      // f-statistic
      // r-squared
      // betas
      // t-stats
      //
      
      // f-statistic
      xoutput[output_index] = REAL(t_sexp)[0];
      
      // r-squared (last value from voxel_lm call: p+2 (stating at 0, so p+1))
      xoutput[output_index + 1 * (sizes[0]*sizes[1]*sizes[2])] = REAL(t_sexp)[p + 1];
      
      // the betas/coefficients:
      for (int k = 2; k < (p + 2); k++) {
        xoutput[output_index + k * (sizes[0]*sizes[1]*sizes[2])] = coefficients[k - 2];
      }
      
      // t-stats
      for(int k = 1; k < p + 1; k++) {
        xoutput[output_index + (k + p + 1) * (sizes[0]*sizes[1]*sizes[2])] = REAL(t_sexp)[k];
      }
  }
   }
 
}	  

}
  Rprintf("\nDone\n");


  for (i=0; i<num_files_left; i++) {
    miclose_volume(hvol_left[i]);
    miclose_volume(hvol_right[i]);
    free(full_buffer_left[i]);
    free(full_buffer_right[i]);
  }

    free(full_buffer_left);
    free(full_buffer_right);
    free(coefficients);
    free(residuals);
    free(effects);
    free(pivot);
    free(work);
    free(qraux);
    free(v);
    free(diag);
    free(se);
    free(t);
    UNPROTECT(3);

  return(output);
}
