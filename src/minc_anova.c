#include "slice_loop_functions.h"


SEXP per_voxel_anova(SEXP filenames, SEXP Sx, SEXP asgn, 
		     SEXP have_mask, SEXP mask) {

  /* generic items for all slice_loop functions */
  unsigned int sizes;
  double **full_buffer;
  double *mask_buffer;
  mihandle_t *hvol, hmask;
  int use_mask;
  int v0, v1, v2;
  int num_files;

  /* anova specific items */
  int n,p, xasgn, maxasgn;
  double *coefficients, *residuals, *effects, *work, *qraux, *v, *diag,
    *se, *t, *comp, *ss;
  int *pivot, *df;
  SEXP Soutput, t_sexp, Sdata;
  double *ouput, *data;

  use_mask = INTEGER(have_mask);
  
  num_files = LENGTH(filenames);
  open_minc_files(filenames, hvol, sizes);
  if (use_mask == 1) {
    get_mask(mask, hmask, mask_buffer, sizes);
  }

  create_slice_buffer(filenames, sizes, full_buffer);

  /* ANOVA specific allocatinos */
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
  PROTECT(Soutput=allocMatrix(REALSXP, (sizes[0] * sizes[1] * sizes[2]), 
			     maxasgn-1));
  output = REAL(Soutput);

  /* allocate the local buffer that will be passed to the function */
  PROTECT(Sdata = allocVector(REALSXP, LENGTH(filenames)));
  data = REAL(Sdata);
  
  /* enter slice loop */
  Rprintf("In slice \n");
  for (v0=0; v0 < sizes[0]; v0++) {
    fill_slice_buffer(filenames, sizes, hvol, buffer, v0);
    if (use_mask == 1) {
      get_mask_slice(hmask, sizes, mask_buffer, v0);
    }
    for (v1=0; v1 < sizes[1]; v1++) {
      for (v2=0; v2 < sizes[2]; v2++) {
	get_indices(v0,v1,v2, sizes, &output_index, &buffer_index);
	
	/* only perform operation if not masked */
	if (use_mask == 0 ||
	    (use_mask == 1 && mask_buffer[buffer_index] > 0.5)) {
	  /* fill data buffer */
	  for (i=0; i< num_files; i++) {
	    data[i] = buffer[i][buffer_index];
	  }

	  /* compute the linear model */
	  t_sexp = voxel_anova(data, Sx, asgn,
			       coefficients, residuals,
			       effects, work, qraux, v, pivot,
			       se, t, comp, ss, df);
	  for (i=0; i < maxasgn-1; i++) {
	    output[output_index + i * (sizes[0]*sizes[1]*sizes[2])] 
	      = REAL(t_sexp)[i];
	  }
	}
	else {
	  for (i=0; i < maxasgn-1; i++) {
	    output[output_index + i * (sizes[0]*sizes[1]*sizes[2])]
	      = 0.0;
	  }
	}
      }
    }
  }
  Rprintf("\nDone\n");
  

  
