#ifndef __MODELLING_FUNCTIONS_H__
#define __MODELLING_FUNCTIONS_H__

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Utils.h>


SEXP voxel_lm(SEXP Sy, SEXP Sx, double *coefficients, 
	      double *residuals, double *effects, 
	      double *work, 
	      double *qraux, double *v, int *pivot, double *se, double *t);

SEXP voxel_lm_file(SEXP Sy, SEXP Sx,int n, int p, double *coefficients, 
	      double *residuals, double *effects, 
	      double *work, 
	      double *qraux, double *v, int *pivot, double *se, double *t);
#endif
