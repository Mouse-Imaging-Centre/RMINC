#include "minc_reader.h"

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
  double x_mean, x_var, y_mean, y_var, x_sum, y_sum, x_n, y_n, s_p;
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

SEXP voxel_lm(SEXP Sy, SEXP Sx, double *coefficients, 
	      double *residuals, double *effects, 
	      double *work, 
	      double *qraux, double *v, int *pivot, double *se, double *t) {

  double tol, rss, resvar;
  double *x, *y, *xoutput;
  int n, p, ny, rank, i, j, rdf, index;
  SEXP new_x, output;


  n = nrows(Sx);
  p = ncols(Sx);

  PROTECT(output=allocVector(REALSXP, p));
  xoutput = REAL(output);

  /* since x is destroyed in dqrls, create a copy here */
  PROTECT(new_x=allocMatrix(REALSXP, n, p));
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

  F77_NAME(dqrls)(x, &n, &p, y, &ny, &tol, coefficients, residuals, effects,
	&rank, pivot, qraux, work);

  //Rprintf("N: %d P: %d NY: %d\n", n,p,ny);
  //Rprintf("Coefficients: ");

  rss = 0;
  rdf = 0;
  for (i=0; i < n; i++) {
    rss += pow(residuals[i], 2);
    //Rprintf("%f\n", residuals[i]);
  }
  rdf = n - p;
  resvar = rss/rdf;
  /*
  for (i=0; i < p; i++) {
    Rprintf("%f  ", coefficients[i]);
  }
  Rprintf("\n");
  Rprintf("RSS: %f   RDF: %d\n", rss, rdf);
  */

  F77_NAME(ch2inv)(x, &n, &p, v);

  //Rprintf("T-stats: ");
  for (i=0; i < p; i++) {
    index = (i * p) + i; /* erm? is this correct for matrix diagonals? */
    se[i] = sqrt(v[index] *resvar);
    //Rprintf("Diag: %f\n", v[index]);
    //Rprintf("SE: %f\n", sqrt(v[index] *resvar));
    xoutput[i] = coefficients[i] / se[i];
    //    Rprintf("%f ", xoutput[i]);
  }
  //Rprintf("\n");

  //Rprintf("R: %f\n", v[0]);


  UNPROTECT(2);
  return(output);
}
  

  
/* minc2_model: run one of a set of modelling function at every voxel
 * filenames: character list of minc2 volumes.
 * y: the variable to model each voxel with. For t-stats and wilcoxon this
 *    is a grouping variable, for correlations it is the correlate. Should
 *    be passed in as a an array of doubles.
 * have_mask: a double of either 0 or 1 depending on whether a mask should
 *            be used.
 * mask: a string containing the mask filename.
 * method: a string containing one of "t-test", "wilcoxon", or "correlation"
 */
SEXP minc2_model(SEXP filenames, SEXP Sx,
			    SEXP have_mask, SEXP mask, SEXP method) {
  int                result;
  mihandle_t         *hvol, hmask;
  char               *method_name;
  int                i, v0, v1, v2, output_index, buffer_index;
  unsigned long      start[3], count[3];
  unsigned long      location[3];
  int                num_files;
  double             *xbuffer, *xoutput, **full_buffer, *xhave_mask, *xn;
  double             *mask_buffer;
  midimhandle_t      dimensions[3];
  unsigned long      sizes[3];
  SEXP               output, buffer, R_fcall, t_sexp;
  /* stuff for linear models only */
  double             *y, *x, *coefficients, *residuals, *effects; 
  double             *diag, *se, *t, *work, *qraux, *v;
  double             tol, rss, resvar;
  int                n, p, ny, rank, j, rdf, index;
  int                *pivot;

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

  /* allocate stuff necessary for fitting linear models */
  if (strcmp(method_name, "lm") == 0) {
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

    PROTECT(t_sexp = allocVector(REALSXP, p));

    /* allocate the output buffer */
    PROTECT(output=allocMatrix(REALSXP, (sizes[0] * sizes[1] * sizes[2]), p));

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
	   || (xhave_mask[0] == 1 && mask_buffer[buffer_index] == 1)) {
	
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
	  else if (strcmp(method_name, "lm") == 0) {
	    t_sexp = voxel_lm(buffer, Sx, coefficients, residuals, effects,
			      work, qraux, v, pivot, se, t);
	    for(i=0; i < p; i++) {
	      xoutput[output_index + i * (sizes[0]*sizes[1]*sizes[2])] 
		      = REAL(t_sexp)[i];
	    }
	  }
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
  else {
    UNPROTECT(2);
  }

  /* return the results */
  return(output);
}
