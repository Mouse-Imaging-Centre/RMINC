#ifndef __MODELLING_FUNCTIONS_H__
#define __MODELLING_FUNCTIONS_H__

//make C functions available from c++
#ifdef __cplusplus 
extern "C" {
#endif

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Utils.h>


SEXP voxel_lm_file(SEXP Sy, SEXP Sx,int n, int p, double *coefficients, 
	      double *residuals, double *effects, 
	      double *work, 
	      double *qraux, double *v, int *pivot, double *se, double *t);

SEXP voxel_lm(SEXP Sy, SEXP Sx,int n,int p,double *coefficients, 
	      double *residuals, double *effects, 
	      double *work, 
	      double *qraux, double *v, int *pivot, double *se, double *t);

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
                 int *df);

#ifdef __cplusplus
}
#endif

#endif
