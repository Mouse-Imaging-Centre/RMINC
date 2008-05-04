#include "slice_loop.h"

/* open_minc_files: given a vector of filenames, open the volumes and
 * returns the dimhandles */
void open_minc_files(SEXP filenames, 
		     midimhandle_t *hvol,
		     unsigned int *sizes) {
  int num_files, i, result;
  midimhandle_t *hvol;

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
  if (result != MI_NOERROR) {
    error("Error getting dimensions\n");
  }
  result = miget_dimension_sizes(dimensions, 3, sizes);
  if (result != MI_NOERROR) {
    error("Error getting dimension sizes\n");
  }

}

/* create_slice_buffer: allocates memory for the slice buffer */
void create_slice_buffer(SEXP filenames,
			 unsigned int *sizes,
			 double **buffer) {
  int num_files, i;

  num_files = LENGTH(filenames);
  buffer = malloc(num_files * sizeof(double));
  for (i=0; i < num_files; i++) {
    buffer[i] = malloc(sizes[1] * sizes[2] * sizeof(double));
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

void free_minc_files(SEXP filenames, midimhandle_t *hvol) {
  int num_files, i;
  
  num_files = LENGTH(filenames);
  for (i=0; i < num_files; i++) {
    miclose_volume(hvol[i]);
    //free
  }
}
