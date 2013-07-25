#include <stdio.h>
//#include <R.h>
//#include <Rinternals.h>
//#include <Rdefines.h>
//#include <R_ext/Utils.h>
#include "modelling_functions.h"

SEXP vertex_lm_loop(SEXP data, SEXP Sx) {

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
  SEXP   output;
  double *xoutput;
  SEXP   buffer;
  double *xbuffer;
  SEXP   t_sexp;
  int    p,n,nVertices,i,j,k;

  xdata = REAL(data);

  // determine numbers of rows, columns, and vertices
  n = nrows(Sx);
  p = ncols(Sx);
  nVertices = nrows(data);
  
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

  // protect voxel_lm output
  PROTECT(t_sexp = allocVector(REALSXP, p+1));;

  // allocate data for output
  PROTECT(output=allocMatrix(REALSXP, nVertices, p+1));
  xoutput=REAL(output);

  // allocate data for the buffer (each vertex for all subjects)
  PROTECT(buffer=allocVector(REALSXP, n));
  xbuffer=REAL(buffer);

  // begin the loop
  Rprintf("Beginning vertex loop: %d %d\n", nVertices, p+1);
  for(i=0; i<nVertices;i++) {
    // fill buffer
    for (j=0; j<n; j++) {
      xbuffer[j] = xdata[i+nVertices*j];
      //Rprintf("B: %d\n", i+nVertices*j);
    }
    t_sexp = voxel_lm(buffer, Sx, coefficients, residuals, effects,
		      work, qraux, v, pivot, se, t);
    for (k=0; k<p+1; k++) {
      //Rprintf("O: %d\n", i+nVertices*k);
      xoutput[i+nVertices*k] = REAL(t_sexp)[k];
    }
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
  UNPROTECT(3);
  
  return(output);
}

SEXP vertex_anova_loop(SEXP data, SEXP Sx,SEXP asgn) {

  double *coefficients;
  double *residuals;
  double *effects;
  int    *pivot;
  int *df;
  double *work;
  double *qraux;
  double *v;
  double *diag;
  double *se;
  double *t;
  double *ss; 
  double *comp;
  double *xdata;
  SEXP   output;
  double *xoutput;
  SEXP   buffer;
  double *xbuffer;
  SEXP   t_sexp;
  int    p,n,nVertices,i,j,k;
  int  maxasgn;
  int *xasgn;
  xdata = REAL(data);

  // determine numbers of rows, columns, and vertices
  n = nrows(Sx);
  p = ncols(Sx);
  nVertices = nrows(data);
  
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
  
  xasgn = INTEGER(asgn);
  maxasgn = 0;
  for (i=0; i < p; i++) {
    if (xasgn[i] > maxasgn) {
      maxasgn = (int) xasgn[i];
    }
  }


  comp = malloc(sizeof(double) * p);
  ss = malloc(sizeof(double) * maxasgn+1);
  df = malloc(sizeof(int) * maxasgn);
  
  
  
  Rprintf("N: %d P: %d\n", n,p);

  // protect voxel_lm output
  PROTECT(t_sexp = allocVector(REALSXP, maxasgn+1));;

  // allocate data for output
  PROTECT(output=allocMatrix(REALSXP, nVertices, maxasgn));
  xoutput=REAL(output);

  // allocate data for the buffer (each vertex for all subjects)
  PROTECT(buffer=allocVector(REALSXP, n));
  xbuffer=REAL(buffer);

  // begin the loop
  Rprintf("Beginning vertex loop: %d %d\n", nVertices, maxasgn+1);
  for(i=0; i<nVertices;i++) {
    // fill buffer
    for (j=0; j<n; j++) {
      xbuffer[j] = xdata[i+nVertices*j];
      //Rprintf("B: %d\n", i+nVertices*j);
    }
    t_sexp = voxel_anova(buffer, Sx,asgn, coefficients, residuals, effects,
		      work, qraux, v, pivot, se, t, comp, ss, df);
    for (k=0; k<maxasgn; k++) {
      //Rprintf("O: %d\n", i+nVertices*k);
      xoutput[i+nVertices*k] = REAL(t_sexp)[k];
    }
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
  free(comp);
  free(ss);

  free(df);
  UNPROTECT(3);
  
  return(output);
}

