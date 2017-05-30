#include "slice_loop.h"
#include "modelling_functions.h"

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


  //Rprintf("ANOVA: N: %d P: %d\n", n,p);
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
  
  //Rprintf("before voxel_lm\n");
  voxel_lm(Sy, Sx,n,p, coefficients, residuals, effects,
	   work, qraux, v, pivot, se, t);
  //Rprintf("after voxel_lm\n");

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

SEXP test_slice_loop(SEXP filenames) {
  misize_t sizes[3];
  double *buffer = NULL;
  mihandle_t *hvol = NULL;
  int i, v0, v1, v2;
  misize_t buffer_index, output_index;
  int num_files;
  SEXP output;
  double *xoutput;


  /*
  Rprintf("allocating memory for volumes\n");
  num_files = LENGTH(filenames);
  hvol = malloc(num_files * sizeof(mihandle_t));
  Rprintf("Before opening files\n");
  open_minc_files(filenames, hvol, sizes);

  Rprintf("Allocating memory for slice buffer\n");
  buffer = malloc(num_files * sizeof(double));
  for (i=0; i < num_files; i++) {
    buffer[i] = malloc(sizes[1] * sizes[2] * sizeof(double));
  }
  */
  Rprintf("Before opening files\n");
  hvol = open_minc_files(filenames, sizes);
  Rprintf("Before allocating slice buffer\n");
  buffer = create_slice_buffer(filenames, sizes);

  Rprintf("protecting output\n");
  PROTECT(output=allocMatrix(REALSXP, (sizes[0] * sizes[1] * sizes[2]), 
			       1));
  xoutput = REAL(output);


  //create_slice_buffer(filenames, sizes, full_buffer);
  Rprintf("In slice \n");
  for (v0=0; v0 < sizes[0]; v0++) {
    Rprintf("Before fill_slice_buffer\n");
    fill_slice_buffer(filenames, sizes, hvol, buffer, v0);
    for (v1=0; v1 < sizes[1]; v1++) {
      for (v2=0; v2 < sizes[2]; v2++) {
	//Rprintf("Before get_indices\n");
	get_indices(v0,v1,v2, sizes, &output_index, &buffer_index);
	//Rprintf("%ul %ul\n", output_index, buffer_index);
	xoutput[output_index] = (double) v0;
      }
    }
  }
  Rprintf("freeing slice buffer\n");
  free_slice_buffer(filenames, sizes, buffer);
  Rprintf("freeing volume handles\n");
  free_minc_files(filenames, hvol);
  Rprintf("finished\n");
  UNPROTECT(1);
  return(output);
}
	 
    

SEXP per_voxel_anova(SEXP filenames, SEXP Sx, SEXP asgn, 
		     SEXP have_mask, SEXP mask, SEXP mask_lower_value, SEXP mask_upper_value) {

  /* generic items for all slice_loop functions */
  misize_t sizes[3];
  double **full_buffer;
  double *mask_buffer;
  mihandle_t *hvol, hmask;
  double *use_mask;
  int v0, v1, v2;
  double xmask_lower_value, xmask_upper_value;
  int num_files, i;
  misize_t buffer_index, output_index;

  /* anova specific items */
  int n,p, maxasgn;
  int *xasgn;
  double *coefficients, *residuals, *effects, *work, *qraux, *v, *diag,
    *se, *t, *comp, *ss;
  int *pivot, *df;
  SEXP Soutput, t_sexp, Sdata;
  double *output, *data;

  use_mask = REAL(have_mask);
  mask_buffer = (double *) calloc(sizes[1]*sizes[2], sizeof(double));
  //Rprintf("Use mask: %f\n", use_mask[0]);
  
  num_files = LENGTH(filenames);
  //Rprintf("Before opening files\n");
  hvol = open_minc_files(filenames, sizes);
  //Rprintf("Before getting mask\n");
  if (use_mask[0] == 1) {
    Rprintf("Getting mask\n");
    //get_mask(mask, hmask, mask_buffer, sizes);
    miopen_volume(CHAR(STRING_ELT(mask, 0)),
		  MI2_OPEN_READ, &hmask);
  }
  
  //Rprintf("Before creating slice buffer\n");
  //full_buffer = malloc(num_files * sizeof(double));
  full_buffer = create_slice_buffer(filenames, sizes);
  //num_files = LENGTH(filenames);
  //full_buffer = malloc(num_files * sizeof(double));
  //for (i=0; i < num_files; i++) {
  //  full_buffer[i] = malloc(sizes[1] * sizes[2] * sizeof(double));
  //}
  /* ANOVA specific allocations */
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
  
  coefficients = (double *) malloc(sizeof(double) * p);
  residuals = (double *) malloc(sizeof(double) * n);
  effects = (double *) malloc(sizeof(double) * n);
  pivot = (int *) malloc(sizeof(int) * p);
  work = (double *) malloc(sizeof(double) * (2*p));
  qraux = (double *) malloc(sizeof(double) * p);
  v = (double *) malloc(sizeof(double) * p * p);
  diag = (double *) malloc(sizeof(double) * p);
  se = (double *) malloc(sizeof(double) * p);
  t = (double *) malloc(sizeof(double) * p);
  
  comp = (double *) malloc(sizeof(double) * p);
  ss = (double *) calloc(maxasgn, sizeof(double));
  df = (int *) malloc(sizeof(int) * maxasgn);
  
  /* get the value at which the mask is to be evaluated */
  xmask_lower_value = REAL(mask_lower_value)[0];
  xmask_upper_value = REAL(mask_upper_value)[0];
  
  Rprintf("N: %d P: %d\n", n,p);
    
  PROTECT(t_sexp = allocVector(REALSXP, maxasgn-1));
  
  //Rprintf("Sizes: %d %d %d\n", sizes[0], sizes[1], sizes[2]);
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
    //Rprintf("Before fill_slice_buffer\n");
    fill_slice_buffer(filenames, sizes, hvol, full_buffer, v0);
    //Rprintf("After fill_slice_buffer\n");
    if (use_mask[0] == 1) {
      get_mask_slice(hmask, sizes, mask_buffer, v0);
    }
    for (v1=0; v1 < sizes[1]; v1++) {
      for (v2=0; v2 < sizes[2]; v2++) {
	//Rprintf("Before get_indices\n");
	get_indices(v0,v1,v2, sizes, &output_index, &buffer_index);
	
	/* only perform operation if not masked */
	if (use_mask[0] == 0 ||
	    (use_mask[0] == 1 && mask_buffer[buffer_index] > (xmask_lower_value - 0.5) && mask_buffer[buffer_index] < (xmask_upper_value + .5))) {
	  /* fill data buffer */
	  for (i=0; i< num_files; i++) {
	    data[i] = full_buffer[i][buffer_index];
	  }

	  /* compute the linear model */
	  //Rprintf("Before voxel_anova\n");
	  t_sexp = voxel_anova(Sdata, Sx, asgn,
			       coefficients, residuals,
			       effects, work, qraux, v, pivot,
			       se, t, comp, ss, df);
	  //Rprintf("Before assigning of output\n");
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
  //Rprintf("freeing slice buffer\n");
  free_slice_buffer(filenames, sizes, full_buffer);
  //Rprintf("freeing volume handles\n");
  free_minc_files(filenames, hvol);
  //Rprintf("finished\n");

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
  Rprintf("\nDone\n");
  return(Soutput);

}
  

  
