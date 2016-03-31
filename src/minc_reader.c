#include "minc_reader.h"

/* compiling:
gcc -shared -fPIC -I/usr/lib/R/include/ -I/projects/mice/share/arch/linux64/include -L/projects/mice/share/arch/linux64/lib -o minc_reader.so minc_reader.c -lminc2 -lhdf5 -lnetcdf -lz -lm
*/

void get_volume_sizes(char **filename, unsigned int *sizes) {
  int result;
  mihandle_t  hvol;
  misize_t tmp_sizes[3];
  midimhandle_t dimensions[3];
   /* open the existing volume */
  result = miopen_volume(filename[0],
			 MI2_OPEN_READ, &hvol);
  if (result != MI_NOERROR) {
    error("Error opening input file: %s.\n", filename[0]);
  }

    /* get the file dimensions and their sizes */
  miget_volume_dimensions( hvol, MI_DIMCLASS_SPATIAL,
			   MI_DIMATTR_ALL, MI_DIMORDER_FILE,
			   3, dimensions);
  result = miget_dimension_sizes( dimensions, 3, tmp_sizes ); 
  //Rprintf("Sizes: %i %i %i\n", tmp_sizes[0], tmp_sizes[1], tmp_sizes[2]);
  sizes[0] = (unsigned int) tmp_sizes[0];
  sizes[1] = (unsigned int) tmp_sizes[1];
  sizes[2] = (unsigned int) tmp_sizes[2];

   // Close the volume, to free handle
   result = miclose_volume(hvol); 
   if (result != MI_NOERROR) {
    error("Error closing file: %s.\n", filename[0]);
  }
  return;
}

SEXP get_vector_from_files(SEXP filenames,  SEXP num_files,  SEXP vec_length,
			  SEXP v1, SEXP v2, SEXP v3) {
  misize_t location[4];
  mihandle_t hvol;
  int result;
  int i, j, vector_length, number_files; 
  SEXP output;
  double *xoutput;

  vector_length = *INTEGER(vec_length);
  number_files = *INTEGER(num_files);

  location[1] = *INTEGER(v1);
  location[2] = *INTEGER(v2);
  location[3] = *INTEGER(v3);

  Rprintf("NFILES: %d NVEC: %d\n", number_files, vector_length);
  PROTECT(output=allocMatrix(REALSXP, number_files, vector_length));
  xoutput = REAL(output);


  for(i=0; i < number_files; i++) {
    /* open the volume */
    result = miopen_volume(CHAR(STRING_ELT(filenames,i)),
			   MI2_OPEN_READ, &hvol);
    if (result != MI_NOERROR) {
      error("Error opening input file: %s.\n", CHAR(STRING_ELT(filenames, i)));
    }

    for(j=0; j < vector_length; j++) {
      location[0] = j;
      result = miget_real_value(hvol, location, 4, &xoutput[i + number_files*j]);

      if (result != MI_NOERROR) {
	error("Error getting voxel from: %s.\n", CHAR(STRING_ELT(filenames, i)));
      }
    }
    miclose_volume(hvol);
  }
  UNPROTECT(1);
  return(output);
}

/* get a voxel from all files, voxel coordinates */
void get_voxel_from_files(char **filenames, int *num_files,
			  int *v1, int *v2, int *v3, double *voxel) {
  misize_t location[3];
  mihandle_t hvol;
  int result;
  int i;

  location[0] = *v1;
  location[1] = *v2;
  location[2] = *v3;

  for(i=0; i < *num_files; i++) {
    /* open the volume */
    result = miopen_volume(filenames[i],
			   MI2_OPEN_READ, &hvol);
    if (result != MI_NOERROR) {
      error("Error opening input file: %s.\n", filenames[i]);
    }

    result = miget_real_value(hvol, location, 3, &voxel[i]);

    if (result != MI_NOERROR) {
      error("Error getting voxel from: %s.\n", filenames[i]);
    }
    miclose_volume(hvol);
  }
}

/* convert voxel coordinates to a vector of world coords */
void convert_voxel_to_world(char **filenames, double *v1, double *v2,
			    double *v3, double *world_coords) {
  int result;
  mihandle_t hvol;
  double voxel_coords[3];

  /* open the volume */
  result = miopen_volume(filenames[0], MI2_OPEN_READ, &hvol);
  if (result != MI_NOERROR) {
    error("Error opening input file: %s.\n", filenames[0]);
  }

  voxel_coords[0] = *v1;
  voxel_coords[1] = *v2;
  voxel_coords[2] = *v3;

  miconvert_voxel_to_world(hvol, voxel_coords, world_coords);
  miclose_volume(hvol);
}

/* convert world coordinates to a vector of voxel coords */
void convert_world_to_voxel(char **filenames, double *v1, double *v2,
			    double *v3, double *voxel_coords) {
  int result;
  mihandle_t hvol;
  double world_coords[3];

  /* open the volume */
  result = miopen_volume(filenames[0], MI2_OPEN_READ, &hvol);
  if (result != MI_NOERROR) {
    error("Error opening input file: %s.\n", filenames[0]);
  }

  world_coords[0] = *v1;
  world_coords[1] = *v2;
  world_coords[2] = *v3;

  miconvert_world_to_voxel(hvol, world_coords, voxel_coords);

  miclose_volume(hvol);
}
  
/* get a voxel from all files, world coordinates */
void get_world_voxel_from_files(char **filenames, int *num_files,
				double *v1, double *v2, double *v3, 
				double *voxel) {
  double location[3], voxel_coord_tmp[3];
  misize_t voxel_coord[3];
  mihandle_t hvol;
  int result;
  int i, j;

  location[0] = *v1;
  location[1] = *v2;
  location[2] = *v3;

  for(i=0; i < *num_files; i++) {
    /* open the volume */
    result = miopen_volume(filenames[i],
			   MI2_OPEN_READ, &hvol);
    if (result != MI_NOERROR) {
      error("Error opening input file: %s.\n", filenames[i]);
    }

    miconvert_world_to_voxel(hvol, location, voxel_coord_tmp);
    
    for (j=0; j < 3; j++) {
      voxel_coord[j] = (unsigned long) voxel_coord_tmp[j] + 0.5;
    }
    result = miget_real_value(hvol, voxel_coord, 3, &voxel[i]);

    if (result != MI_NOERROR) {
      error("Error getting voxel from: %s.\n", filenames[i]);
    }
    miclose_volume(hvol);
  }
}

SEXP minc_history_size(SEXP filename){
  size_t hist_size;
  mihandle_t hvol;
  int result;
  const char *filepath = CHAR(asChar(filename));
  
  result = miopen_volume(filepath,
                         MI2_OPEN_READ, &hvol);
  
  if (result != MI_NOERROR) {
    error("Error opening input file: %s.\n", filepath);
  }
  
  miget_attr_length(hvol, "", "history", &hist_size);

  miclose_volume(hvol);
  
  SEXP output = ScalarInteger((int) hist_size);
  return(output);
}

SEXP get_minc_history(SEXP filename) {
  int                result;
  mihandle_t         hvol;
  int int_size   =   asInteger(minc_history_size(filename));
  char  *history = malloc(int_size);  
  const char *filepath = CHAR(asChar(filename));

  result = miopen_volume(filepath,
                         MI2_OPEN_READ, &hvol);

  if (result != MI_NOERROR) {
    error("Error opening input file: %s.\n", filepath);
  }
  
  miget_attr_values(hvol, MI_TYPE_STRING, "", "history", int_size, history);
  miclose_volume(hvol);
  
  SEXP output = PROTECT(mkString(history));
  free(history);
  UNPROTECT(1);
  
  return(output);
}

SEXP minc_overwrite_history(SEXP filename, SEXP history, SEXP hist_size){
  const char *history_line = CHAR(asChar(history));
  const char *filepath = CHAR(asChar(filename));
  int history_size = asInteger(hist_size);
  mihandle_t  hvol;
  int read_result;
  int hist_edit_result;
  
  read_result = miopen_volume(filepath,
                              MI2_OPEN_RDWR, &hvol);
  
  if (read_result != MI_NOERROR) {
    error("Error opening input file: %s.\n", filepath);
  }
  
  hist_edit_result = miadd_history_attr(hvol, history_size, history_line);
  if(hist_edit_result != MI_NOERROR){
    error("Error editing history for file: %s \n", filepath);
  }
  
  miclose_volume(hvol);
  
  return(R_NilValue);
}
  

/* get a real value hyperslab from file */
void get_hyperslab(char **filename, int *start, int *count, double *slab) {
  int                result;
  mihandle_t         hvol;
  int                i;
  unsigned long      tmp_start[3];
  unsigned long      tmp_count[3];

  /* open the volume */
  result = miopen_volume(filename[0],
			 MI2_OPEN_READ, &hvol);
  if (result != MI_NOERROR) {
    error("Error opening input file: %s.\n", filename[0]);
  }

  for (i=0; i < 3; i++) {
    tmp_start[i] = (unsigned long) start[i];
    tmp_count[i] = (unsigned long) count[i];
  }

  /* get the hyperslab */
  Rprintf("Start: %i %i %i\n", start[0], start[1], start[2]);
  Rprintf("Count: %i %i %i\n", count[0], count[1], count[2]);
  if (miget_real_value_hyperslab(hvol, 
				 MI_TYPE_DOUBLE, 
				 (misize_t *) tmp_start, 
				 (misize_t *) tmp_count, 
				 slab)
      < 0) {
    error("Could not get hyperslab.\n");
  }
   // Close the volume, to free handle
   result = miclose_volume(hvol); 
   if (result != MI_NOERROR) {
    error("Error closing file: %s.\n", filename[0]);
  }

  return;
}

/* get a real value hyperslab from file */
SEXP get_hyperslab2( SEXP filename,  SEXP start,  SEXP count, SEXP slab) {
  int                result;
  mihandle_t         hvol;
  int                i;
  misize_t      tmp_start[3];
  misize_t      tmp_count[3];

  /*
  char **c_filename;
  int *c_start;
  int *c_count;
  double *c_slab;
  */
  /* open the volume */
  
  Rprintf("Crap %s\n", CHAR(STRING_ELT(filename, 0)));
  result = miopen_volume(CHAR(STRING_ELT(filename,0)),
			 MI2_OPEN_READ, &hvol);
  if (result != MI_NOERROR) {
    error("Error opening input file: %s.\n", CHAR(STRING_ELT(filename,0)));
  }

  for (i=0; i < 3; i++) {
    tmp_start[i] = (unsigned long) INTEGER(start)[i];
    tmp_count[i] = (unsigned long) INTEGER(count)[i];
  }

  /* get the hyperslab */
  Rprintf("Start: %i %i %i\n", INTEGER(start)[0], INTEGER(start)[1], INTEGER(start)[2]);
  Rprintf("Count: %i %i %i\n", INTEGER(count)[0], INTEGER(count)[1], INTEGER(count)[2]);
  if (miget_real_value_hyperslab(hvol, 
				 MI_TYPE_DOUBLE, 
				 tmp_start, 
				 tmp_count, 
				 REAL(slab))
      < 0) {
    error("Could not get hyperslab.\n");
  }
  return(slab);
}

  

/* minc2_apply: evaluate a function at every voxel of a series of files
 * filenames: character array of filenames. Have to have identical sampling.
 * fn: string representing a function call to be evaluated. The variable "x"
 *     will be a vector of length number_volumes in the same order as the 
 *     filenames array.
 * have_mask: 0 if there is no mask, 1 if there is
 * mask: filename containing the mask
 * mask_value: value in the mask where function is to be evaluated
 * rho: the R environment.
 */
     
SEXP minc2_apply(SEXP filenames, SEXP fn, SEXP have_mask, 
		 SEXP mask, SEXP mask_value, SEXP rho) {
  int                result;
  mihandle_t         *hvol, hmask;
  int                i, v0, v1, v2, output_index, buffer_index;
  unsigned long     start[3], count[3];
  //unsigned long      location[3];
  int                num_files;
  double             *xbuffer, *xoutput, **full_buffer;
  double             *xhave_mask, *xmask_value;
  double             *mask_buffer;
  midimhandle_t      dimensions[3];
  misize_t            sizes[3];
  SEXP               output, buffer;
  //SEXP               R_fcall;
  

  /* allocate memory for volume handles */
  num_files = LENGTH(filenames);
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
  
  /* get the value inside that mask */
  xmask_value = REAL(mask_value);

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
				     (misize_t *) start, 
				     (misize_t *) count, 
				     full_buffer[i]) )
	error("Error opening buffer.\n");
    }
    /* get mask - if desired */
    if (xhave_mask[0] == 1) {
      if (miget_real_value_hyperslab(hmask, 
				     MI_TYPE_DOUBLE, 
				     (misize_t *) start, 
				     (misize_t *) count, 
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
	       mask_buffer[buffer_index] > xmask_value[0] -0.5 &&
	       mask_buffer[buffer_index] < xmask_value[0] + 0.5)) {
	
	  for (i=0; i < num_files; i++) {
// 	    location[0] = v0;
// 	    location[1] = v1;
// 	    location[2] = v2;
	    //SET_VECTOR_ELT(buffer, i, full_buffer[i][index]);
	    //result = miget_real_value(hvol[i], location, 3, &xbuffer[i]);
	    xbuffer[i] = full_buffer[i][buffer_index];
	    
	    //Rprintf("V%i: %f\n", i, full_buffer[i][index]);

	  }
	  /* install the variable "x" into environment */
	  defineVar(install("x"), buffer, rho);
	  //SETCADDR(R_fcall, buffer);
	  //SET_VECTOR_ELT(output, index, eval(R_fcall, rho));
	  //SET_VECTOR_ELT(output, index, test);
	  /* evaluate the function */
	  xoutput[output_index] = REAL(eval(fn, rho))[0]; 
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
  UNPROTECT(2);

  /* return the results */
  return(output);
}

/* writes a hyperslab to the output filename, creating the output voluem
   to be like the second filename passed in */
void write_minc2_volume(char **output, char **like_filename,
			int *start, int *count, double *max_range,
			double *min_range, double *slab) {
  mihandle_t hvol_like, hvol_new;
  midimhandle_t *dimensions_like, *dimensions_new;
  unsigned long tmp_count[3];
  unsigned long tmp_start[3];
  int i;

  /* allocate the dimension handles */
  dimensions_like = malloc(sizeof(midimhandle_t) * 3);
  dimensions_new = malloc(sizeof(midimhandle_t) * 3);

  /* read the like volume */
  if (miopen_volume(like_filename[0], MI2_OPEN_READ, &hvol_like) < 0 ) {
    error("Error opening volume: %s\n", like_filename[0]);
  }
  
  /* get dimensions */
  if (miget_volume_dimensions( hvol_like, MI_DIMCLASS_SPATIAL, MI_DIMATTR_ALL,
			       MI_DIMORDER_FILE, 3, dimensions_like ) < 0 ) {
    error("Error getting volume dimensions\n");
  }

  /* copy the dimensions to the new file */
  for (i=0; i < 3; i++) {
    if (micopy_dimension(dimensions_like[i], &dimensions_new[i]) < 0) {
      error ("Error copying dimension %d\n", i);
    }
  }

  /* create the new volume */
  if ( micreate_volume(output[0], 3, dimensions_new, MI_TYPE_USHORT,
		       MI_CLASS_REAL, NULL, &hvol_new) < 0 ) {
    error("Error creating volume %s\n", output[0]);
  }
  if (micreate_volume_image(hvol_new) < 0) {
    error("Error creating volume image\n");
  }
  
  /* set valid and real range */
  miset_volume_valid_range(hvol_new, 65535, 0);
  miset_volume_range(hvol_new, max_range[0], min_range[0]);

  Rprintf("Range: %f %f\n", max_range[0], min_range[0]);

  /* write the buffer */
  for (i=0; i < 3; i++) {
    tmp_start[i] = (unsigned long) start[i];
    tmp_count[i] = (unsigned long) count[i];
  }
  if (miset_real_value_hyperslab(hvol_new, MI_TYPE_DOUBLE, 
				 (misize_t *) tmp_start, 
				 (misize_t *) tmp_count,
				 slab) < 0) {
    error("Error writing buffer to volume\n");
  }
  if (miclose_volume(hvol_like) < 0) {
    error("Error closing like volume\n");
  }
  if (miclose_volume(hvol_new) < 0) {
    error("Error closing new volume\n");
  }

  free(dimensions_new);
  free(dimensions_like);
  return;
}

