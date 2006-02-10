#include "minc_reader.h"


/* compute a  t test given a voxel and grouping */
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

SEXP voxel_correlation(SEXP Sx, SEXP Sy, SEXP Sn) {
  double *x, *y, *n, *r;
  double sum_x, sum_y, sum_xy, sum_x2, sum_y2, numerator, denominator;
  int i;
  SEXP output;

  PROTECT(output=allocVector(REALSXP, 1));
  
  r = REAL(output);
  x = REAL(Sx);
  y = REAL(Sy);
  n = REAL(Sn);

  /* compute sums for x and y */
  sum_x = 0;
  sum_y = 0;
  sum_xy = 0;
  sum_y2 = 0;
  sum_x2 = 0;
  for (i=0; i < *n; i++) {
    sum_x += x[i];
    sum_y += y[i];
    sum_xy += x[i] * y[i];
    sum_x2 += pow(x[i], 2);
    sum_y2 += pow(y[i], 2);
  }

  numerator = *n * sum_xy - sum_x * sum_y;
  denominator = sqrt( (*n * sum_x2 - pow(sum_x, 2)) * 
		      (*n * sum_y2 - pow(sum_y, 2)) );
  *r = numerator / denominator;
  UNPROTECT(1);
  return(output);
}

  
     
SEXP minc2_group_comparison(SEXP filenames, SEXP groupings, SEXP na, SEXP nb,
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
  SEXP               output, buffer, R_fcall, n;
  

  /* determine the number of files */
  PROTECT(n=allocVector(REALSXP, 1));
  xn = REAL(n);
  *xn = LENGTH(filenames);

  num_files = (int) *xn;

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

  /* allocate the output buffer */
  PROTECT(output=allocVector(REALSXP, (sizes[0] * sizes[1] * sizes[2])));
  xoutput = REAL(output);

  /* allocate the local buffer that will be passed to the function */
  PROTECT(buffer=allocVector(REALSXP, num_files));
  xbuffer = REAL(buffer); 

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
	    xoutput[output_index] = REAL(t_test(buffer, groupings, n))[0]; 
	  }
	  else if (strcmp(method_name, "wilcoxon") == 0) {
	    xoutput[output_index] = 
	      REAL(wilcoxon_rank_test(buffer, groupings, na, nb))[0];
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
  free(full_buffer);
  UNPROTECT(3);

  /* return the results */
  return(output);
}
