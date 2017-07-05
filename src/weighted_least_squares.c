#include <stdio.h>
#include "modelling_functions.h"
#include "R_ext/Lapack.h"
#include "R_ext/Applic.h"

SEXP voxel_wlm(SEXP Sy, SEXP Sx, SEXP ws, int n,int p,double *coefficients, 
              double *residuals, double *effects, 
              double *work, 
              double *qraux, double *v, int *pivot, double *se, double *t) {
  
  double tol, rss, resvar, mss, mean_fitted, sum_fitted, logLik;
  double *x, *y, *xoutput, *xws;
  int ny,ny1, rank, i, j, rdf, index;
  SEXP new_x, output;
  int nprot = 0;
  
  
  ny = ncols(Sy);
  ny1 = nrows(Sy);
  
  // the output will contain:
  // 
  // f-statistic
  // p * t-statistic
  // r-squared
  // logLik
  //
  // which is a total of p+3 values
  PROTECT(output=allocVector(REALSXP, p+3)); nprot++;
  xoutput = REAL(output);
  
  double total_log_weight = 0;
  xws = (double *) malloc(n * sizeof(double));
  for(i = 0; i < n; ++i){
    total_log_weight += log(REAL(ws)[i]);
    xws[i] = sqrt(REAL(ws)[i]);
    //Rprintf("%f\n", xws[i]);
  }
  /* since x (the model matrix, input variable Sx) is destroyed in dqrls, 
  create a copy here */
  PROTECT(new_x=allocMatrix(REALSXP, n, p)); nprot++;
  for (i=0; i < n; i++) {
    for (j=0; j < p; j++) {
      REAL(new_x)[i+j*n] = REAL(Sx)[i+j*n] * xws[i]; //multiply in weights
      //Rprintf("new_x %f\n", REAL(new_x)[i+j*n]);
    }
  }
  
  x = REAL(new_x);
  y = (double *) malloc(n * sizeof(double));
  for(i=0; i < n; ++i)
    y[i] = REAL(Sy)[i] * xws[i];
  
  //Rprintf("coly %d rowy %d\n", ny,ny1);
  rank = 1;
  tol = 1e-07;
  for(int j = 0; j < p; ++j)
    pivot[j] = j;
  
  // compute the least squares solution:
  F77_CALL(dqrls)(x, &n, &p, y, &ny, &tol, coefficients, residuals, effects,
           &rank, pivot, qraux, work);
  
  double *piv_coefficients = (double *) malloc(p * sizeof(double));
  memcpy(piv_coefficients, coefficients, p * sizeof(double));
  for(j = 0; j < p; ++j)
    coefficients[pivot[j]] = piv_coefficients[j];
  free(piv_coefficients);
  
  // Calculate the f-statistic first
  rss = 0; // residual sum of squares
  double wrss = 0; //weighted residual sum of squares
  rdf = 0; // residual degrees of freedom
  mss = 0; // (fitted - mean fitted sum) of squares
  sum_fitted = 0;
  logLik = 0;
  for (i=0; i < n; i++) {
    residuals[i] /= xws[i];
    double res2 = pow(residuals[i],2);
    rss += res2;
    wrss += REAL(ws)[i] * res2; 
    sum_fitted += (y[i] - residuals[i]);
  }
  
  logLik = .5 * (total_log_weight - n * (log(2*M_PI) + 1 - log(n) + log(wrss)));
  
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
  //Rprintf("in voxel_lm, F: %f\n", xoutput[0]);
  
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
  
  // the r-squared:
  xoutput[p+1] = mss / (mss + rss);
  xoutput[p+2] = logLik;
  free(xws);
  free(y);
  UNPROTECT(nprot);
  return(output);
}

SEXP vertex_wlm_loop(SEXP data_left, SEXP data_right,SEXP mmatrix, SEXP ws)  {
  
  double *coefficients;
  double *residuals;
  double *effects;
  int    *pivot;
  double *work;
  double *qraux;
  double *v;
  double *diag;
  double *se;
  double *t;
  double *xdata;
  double *ydata;
  double *pMmatrix;
  SEXP   output;
  double *xoutput;
  SEXP   buffer;
  double *xbuffer;
  SEXP   buffer1;
  double *ybuffer;  
  SEXP   t_sexp;
  int    p,n,nVertices,i,j,k,mmatrix_rows,mmatrix_cols;
  
  xdata = REAL(data_left);
  if(!isLogical(data_right)) 
    ydata = REAL(data_right);
  //Rprintf("Starting\n");
  // Case 1: There is no static part
  if(isLogical(mmatrix)) 
    // For now, maximum allowed dynamic parts is 1 so set p = 2 (1 for intercept)
    p = 2;
  // Case 2: There is a static part
  else {
    if(isLogical(data_right)) {
      //Rprintf("in isLogical(right)\n");
      pMmatrix = REAL(mmatrix);
      mmatrix_cols = ncols(mmatrix);
      mmatrix_rows = nrows(mmatrix);
      p =  mmatrix_cols;
    }
    else {
      pMmatrix = REAL(mmatrix);
      mmatrix_cols = ncols(mmatrix);
      mmatrix_rows = nrows(mmatrix);
      p = mmatrix_cols + 1;
      //Rprintf("mmatrix cols: %d mmatrix rows: %d\n", mmatrix_cols,mmatrix_rows );
    }
  }
  
  n = ncols(data_left);
  nVertices = nrows(data_left);
  
  
  // allocate memory for lm
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

  // protect voxel_lm output, currently returns:
  // 
  // f-statistic
  // p * t-statistic
  // r-squared
  // 
  //Rprintf("size of t_sexp: %d\n", p+2);
  PROTECT(t_sexp = allocVector(REALSXP, p + 3));
  //Rprintf("output 0: %f\n", REAL(t_sexp)[0]);
  // allocate data for output, will contain:
  // 
  // f-statistic
  // r-squared
  // betas
  // t-stats
  // 
  PROTECT(output=allocMatrix(REALSXP, nVertices, 2*p + 3));
  xoutput=REAL(output);
  
  // allocate data for the buffer (each vertex for all subjects)
  PROTECT(buffer=allocVector(REALSXP, n));
  xbuffer=REAL(buffer);
  
  PROTECT(buffer1=allocVector(REALSXP, n*p));
  ybuffer=REAL(buffer1);
  
  // begin the loop
  Rprintf("Beginning vertex loop: %d %d\n", nVertices, p+1);
  for(i=0; i<nVertices;i++) {
    // fill buffer
    for (j=0; j<n; j++) {
      xbuffer[j] = xdata[i+nVertices*j];
      //Rprintf("B: %d\n", i+nVertices*j);
      //Rprintf("data: %f\n", xbuffer[j]);
    }
    // If no mmatrix have to build it ourselves 
    if(isLogical(mmatrix))  { 
      // fill y buffer
      // Intercept
      for (j=0; j < n; j++) {
        ybuffer[j] = 1.0;
        //Rprintf("ybuffer %f index %d\n", ybuffer[j] ,j);
      }    
      // Current Vertex Data
      for (j=0; j<n; j++) {
        ybuffer[j+n] = ydata[i+nVertices*j];
        //Rprintf("ybuffer %f index %d\n", ybuffer[j+n],j+n);
      }    
    }
    else {
      //Rprintf("filling with static part\n");
      // Fill with static part
      for (j=0; j < mmatrix_cols*mmatrix_rows; j++) {
        ybuffer[j] = pMmatrix[j];
        //Rprintf("mmatrix %f index %d\n", pMmatrix[j],j);
      }  
      // Fill with dynamic part  (if exists)
      if(!isLogical(data_right)) {
        //Rprintf("!isLogical(right)\n");
        for (j=0; j < n ; j++) {
          ybuffer[j+mmatrix_cols*mmatrix_rows] = ydata[i+nVertices*j];
          //Rprintf("mmatrix %f index %d\n", ydata[i+nVertices*j],j+mmatrix_cols*mmatrix_rows);
        }   
      }
    }
    t_sexp = voxel_wlm(buffer,buffer1,ws,n,p,coefficients, residuals, effects,
                      work, qraux, v, pivot, se, t);
    // f-statistic
    xoutput[i] = REAL(t_sexp)[0];
    // r-squared (last value from voxel_lm call: p+2 (stating at 0, so p+1))
    xoutput[i + 1 * nVertices] = REAL(t_sexp)[p+1];
    // the betas/coefficients:
    for (int k = 2; k < (p + 2); k++) {
      xoutput[i + k * nVertices] = coefficients[k - 2];
    }
    // t-stats
    for(int k = 1; k < p + 1; k++) {
      xoutput[i + (k + p + 1) * nVertices] = REAL(t_sexp)[k];
    }
    
    //logLik
    xoutput[i + (2 * p + 2) * nVertices] = REAL(t_sexp)[p+2];
  }
  
  Rprintf("Done with vertex loop\n");
  
  // free memory
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
  UNPROTECT(4);
  
  return(output);
}