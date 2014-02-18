#include <stdio.h>
//#include <R.h>
//#include <Rinternals.h>
//#include <Rdefines.h>
//#include <R_ext/Utils.h>
#include "modelling_functions.h"

SEXP vertex_lm_loop(SEXP data_left, SEXP data_right,SEXP mmatrix)  {

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

  // Case 1: There is no static part
  if(isLogical(mmatrix)) 
	// For now, maximum allowed dynamic parts is 1 so set p = 2 (1 for intercept)
	p = 2;
  // Case 2: There is a static part
  else {
 	if(isLogical(data_right)) 
{
         pMmatrix = REAL(mmatrix);
  	 mmatrix_cols = ncols(mmatrix);
  	 mmatrix_rows = nrows(mmatrix);
  	  p =  mmatrix_cols ;}
else {
        pMmatrix = REAL(mmatrix);
  	mmatrix_cols = ncols(mmatrix);
  	mmatrix_rows = nrows(mmatrix);
        p = mmatrix_cols + 1;
        Rprintf("mmatrix cols: %d mmatrix rows: %d\n", mmatrix_cols,mmatrix_rows ); }
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
  PROTECT(t_sexp = allocVector(REALSXP, p + 2));;

  // allocate data for output, will contain:
  // 
  // f-statistic
  // r-squared
  // betas
  // t-stats
  // 
  PROTECT(output=allocMatrix(REALSXP, nVertices, 2*p + 2));
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
		// Fill with static part
	    	for (j=0; j < mmatrix_cols*mmatrix_rows; j++) {
			ybuffer[j] = pMmatrix[j];
		        //Rprintf("mmatrix %f index %d\n", pMmatrix[j],j);
		    }  
		// Fill with dynamic part  (if exists)
 		if(!isLogical(data_right)) {
		    	for (j=0; j < n ; j++) {
				ybuffer[j+mmatrix_cols*mmatrix_rows] = ydata[i+nVertices*j];
				//Rprintf("mmatrix %f index %d\n", ydata[i+nVertices*j],j+mmatrix_cols*mmatrix_rows);
			    }   
		}
	    }

    t_sexp = voxel_lm(buffer,buffer1,n,p,coefficients, residuals, effects,
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
  maxasgn++;	

  comp = malloc(sizeof(double) * p);
  ss = malloc(sizeof(double) * maxasgn);
  df = malloc(sizeof(int) * maxasgn);
  
  Rprintf("N: %d P: %d\n", n,p);

  // protect output
  PROTECT(t_sexp = allocVector(REALSXP, maxasgn));

  // allocate data for output
  PROTECT(output=allocMatrix(REALSXP, nVertices, maxasgn-1));
  xoutput=REAL(output);

  // allocate data for the buffer (each vertex for all subjects)
  PROTECT(buffer=allocVector(REALSXP, n));
  xbuffer=REAL(buffer);

  // begin the loop
  Rprintf("Beginning vertex loop: %d %d\n", nVertices, maxasgn);
  for(i=0; i<nVertices;i++) {
    // fill buffer
    for (j=0; j<n; j++) {
      xbuffer[j] = xdata[i+nVertices*j];
    }
    t_sexp = voxel_anova(buffer, Sx,asgn, coefficients, residuals, effects,
		      work, qraux, v, pivot, se, t, comp, ss, df);
    for (k=0; k<maxasgn-1; k++) {
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

