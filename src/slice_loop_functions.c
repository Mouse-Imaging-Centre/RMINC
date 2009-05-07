#include "slice_loop.h"

/* open_minc_files: given a vector of filenames, open the volumes and
 * returns the dimhandles */
mihandle_t* open_minc_files(SEXP filenames, 
			       unsigned int *sizes) {
  int num_files, i, result;
  midimhandle_t dimensions[3];
  mihandle_t *hvol;

  num_files = LENGTH(filenames);
  hvol = malloc(num_files * sizeof(mihandle_t));

  /* open all the files */
  for (i=0; i < num_files; i++) {
    result = miopen_volume(CHAR(STRING_ELT(filenames, i)),
			   MI2_OPEN_READ, &hvol[i]);
    if (result != MI_NOERROR) {
      error("Error opening input file: %s.\n", CHAR(STRING_ELT(filenames,i)));
    }
  }

  /* get dimension size - assume it's the same for all volumes */
  result = miget_volume_dimensions( hvol[0], MI_DIMCLASS_SPATIAL,
				    MI_DIMATTR_ALL, MI_DIMORDER_FILE,
				    3, dimensions );
  if (result == MI_ERROR) {
    error("Error getting dimensions\n");
  }
  result = miget_dimension_sizes(dimensions, 3, sizes);
  Rprintf("Sizes: %d %d %d\n", sizes[0], sizes[1], sizes[2]);
  if (result == MI_ERROR) {
    error("Error getting dimension sizes\n");
  }
  return(hvol);

}

/* get_mask: opens the mask volume and allocates memory for the 
 * mask buffer */
void get_mask(SEXP filename, mihandle_t hmask, double *mask_buffer,
	      unsigned int *sizes) {
  int result;
  
  result = miopen_volume(CHAR(STRING_ELT(filename, 0)),
			 MI2_OPEN_READ, &hmask);
  if (result != MI_NOERROR) {
    error("Error opening mask: %s.\n", CHAR(STRING_ELT(filename, 0)));
  }
  mask_buffer = malloc(sizes[1] * sizes[2] * sizeof(double));
}

/* create_slice_buffer: allocates memory for the slice buffer */
double** create_slice_buffer(SEXP filenames,
			     unsigned int *sizes) {

  int num_files, i;
  double **buffer;
  
  //Rprintf("before: %p\n", buffer[0]);

  num_files = LENGTH(filenames);
  buffer = malloc(num_files * sizeof(double));
  for (i=0; i < num_files; i++) {
    buffer[i] = malloc(sizes[1] * sizes[2] * sizeof(double));
  }
  //Rprintf("after: %p %d %d\n", buffer[0], sizes[1], sizes[2]);
  return(buffer);
}

void get_indices(int v0, int v1, int v2, unsigned int *sizes,
		 unsigned long *output_index, 
		 unsigned long *buffer_index) {
  *output_index = (unsigned long) v0*sizes[1]*sizes[2]+v1*sizes[2]+v2;
  *buffer_index = (unsigned long) sizes[2] * v1 + v2;
}


void fill_slice_buffer(SEXP filenames, 
		       unsigned int *sizes,
		       mihandle_t *hvol,
		       double **buffer,
		       int slice_number) {
  int num_files, result, i;
  unsigned long start[3];
  unsigned long count[3];

  //Rprintf("Sizes in fill_slice: %d %d %d\n",
  //sizes[0], sizes[1], sizes[2]);
  Rprintf("%d ", slice_number);

  start[0] = slice_number; start[1] = 0; start[2] = 0;
  count[0] = (unsigned long) 1; 
  count[1] = (unsigned long) sizes[1]; 
  count[2] = (unsigned long) sizes[2];
  //count[1] = 1;
  //count[2] = 1;

  num_files = LENGTH(filenames);
  for (i=0; i < num_files; i++) {
    //Rprintf("buffer: %p\n", buffer[0]);
    //Rprintf("fill_slice_buffer: f %d, start %lu %lu %lu, count %lu %lu %lu\n",
    //i, start[0], start[1], start[2], 
    //count[0], count[1], count[2]);
    result = miget_real_value_hyperslab(hvol[i], 
					MI_TYPE_DOUBLE,
					(unsigned long *)start, 
					(unsigned long *)count, 
					buffer[i]);
    //Rprintf("hs results: %d\n", result);
    if (result != MI_NOERROR) {
      error("Error getting data from slice %d.\n", slice_number);
    }
  }
}

void get_mask_slice(mihandle_t hmask,
		    unsigned int *sizes,
		    double *mask_buffer,
		    int slice_number) {
  int result;
  unsigned long start[3];
  unsigned long count[3];

  start[0] = slice_number; start[1] = 0; start[2] = 0;
  count[0] = 1; count[1] = sizes[1]; count[2] = sizes[2];

  result = miget_real_value_hyperslab(hmask, MI_TYPE_DOUBLE,
				      start, count, mask_buffer);
  if (result != MI_NOERROR) {
    error("Error getting mask slice %d.\n", slice_number);
  }
}

void free_slice_buffer(SEXP filenames,
		       unsigned int *sizes,
		       double **buffer) {
  int num_files, i;
  
  num_files = LENGTH(filenames);
  for (i=0; i < num_files; i++) {
    free(buffer[i]);
  }
}

void free_minc_files(SEXP filenames, mihandle_t *hvol) {
  int num_files, i;
  
  num_files = LENGTH(filenames);
  for (i=0; i < num_files; i++) {
    miclose_volume(hvol[i]);
    //free
  }
}
