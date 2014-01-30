#include "minc2_io.h"



/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * Purpose: Given a set of 3-D coordinates in x,y,z order, traverse
 *          through all frames of all specified files, and return a  
 *          matrix containing all of the real values at the given voxel
 *          coordinate.
 *
 *          The dimensions of the returned matrix is:
 *               (number_of_files * number_of_frames) 
 *
 *          Note that passing a single 3-D volume will return a 1x1 matrix.
*/
SEXP read_voxel_from_files(SEXP filenames,  SEXP voxCoords,  SEXP noFiles,  SEXP noFrames) {
	mihandle_t		minc_volume;
	int				result;
	int				n_dimensions;
	int 			i, ndx;
	int				no_files; 
	int				no_frames;
	int				no_rows;
	int				no_cols;
	unsigned long	hSlab_start[MI2_MAX_VAR_DIMS];
	unsigned long	hSlab_count[MI2_MAX_VAR_DIMS];
	int				hSlab_buffer_size;
	static char *dimorder3d[] = { "zspace","yspace","xspace" };
	static char *dimorder4d[] = { "time", "zspace","yspace","xspace" };
	
	
	SEXP			output;
	int				output_ndx;
	double			*hSlab_buffer;


	// start ...
	if ( R_DEBUG_mincIO ) Rprintf("read_voxel_from_files: start ...\n");


	// init
	no_files = INTEGER(noFiles)[0];
	no_frames = INTEGER(noFrames)[0];


	/*set the xyz coords of the sampled point
	... hSlab_start[0] will hold the frame index */
	hSlab_start[0] = 0;
	hSlab_start[1] = INTEGER(voxCoords)[0];
	hSlab_start[2] = INTEGER(voxCoords)[1];
	hSlab_start[3] = INTEGER(voxCoords)[2];


	// allocate receiving matrix (3-D volumes always only return 1 column)
	no_rows = no_files;
	no_cols = (no_frames == 0) ? 1 : no_frames;
	PROTECT(output=allocMatrix(REALSXP, no_rows, no_cols));


	// allocate the hyper-slab output buffer to the heap (big enough to handle ONE hyper-slab)
	hSlab_buffer_size =  no_cols;
	if ( R_DEBUG_mincIO ) Rprintf("Allocating hyperslab buffer using R_alloc.  Size is %d doubles\n", hSlab_buffer_size);
	hSlab_buffer = (double *) R_alloc(hSlab_buffer_size, sizeof(double));
	if ( hSlab_buffer == NULL) {
		error("Error allocating hyperslab buffer.  Size: %d doubles.\n", hSlab_buffer_size);
	}


	// loop over all volumes and extract a single value for each frame in the volume
	for( i=0; i < no_files; ++i ) {

		/* open the volume */
		result = miopen_volume(CHAR(STRING_ELT(filenames,i)), MI2_OPEN_READ, &minc_volume);
		//
		if (result != MI_NOERROR) {
			error("Error opening input file: %s.\n", CHAR(STRING_ELT(filenames, i)));
		}


		/* set the apparent order to something conventional */
		//	... first need to get the number of dimensions
		if ( miget_volume_dimension_count(minc_volume, MI_DIMCLASS_ANY, MI_DIMATTR_ALL, &n_dimensions) != MI_NOERROR ){
			error("Error returned from miget_volume_dimension_count.\n");
		}
		/* ... now set the order */
		if ( R_DEBUG_mincIO ) Rprintf("read_voxel_from_files: Setting the apparent order for %d dimensions ... ", n_dimensions);
		if ( n_dimensions == 3 ) {
			result = miset_apparent_dimension_order_by_name(minc_volume, 3, dimorder3d);
		} else if ( n_dimensions == 4 ) {
			result = miset_apparent_dimension_order_by_name(minc_volume, 4, dimorder4d);
		} else {
			error("Error file %s has %d dimensions and we can only deal with 3 or 4.\n", CHAR(STRING_ELT(filenames,i)), n_dimensions);
		}

		if ( result != MI_NOERROR ) { 
			error("Error returned from miset_apparent_dimension_order_by_name while setting apparent order for %d dimensions.\n", n_dimensions); 
		}
		if ( R_DEBUG_mincIO ) Rprintf("Done.\n");



		if ( no_frames > 0 ) {

			// read a hyperslab across all frames (i.e. over time)
			hSlab_count[0] = no_frames;
			hSlab_count[1] = hSlab_count[2] = hSlab_count[3] = 1;
			if ( R_DEBUG_mincIO ) Rprintf("hSlab_count [0..3] = %d, %d, %d, %d\n", 
							hSlab_count[0], hSlab_count[1], hSlab_count[2], hSlab_count[3]);
			
			result = miget_real_value_hyperslab(minc_volume, MI_TYPE_DOUBLE, 
																hSlab_start,
																hSlab_count,
																hSlab_buffer);
			if ( result != MI_NOERROR ) {
				error("Error in miget_real_value_hyperslab: %s.\n", CHAR(STRING_ELT(filenames, i)));
			}
			
			// move values from hyper-slab buffer to output buffer
			// ... we're doing this because R expects data in column-major order, 
			// ... and we're writing the frame values as rows
			for ( ndx=0; ndx < no_frames; ++ndx) {
				output_ndx = (no_files * ndx) +i;
				REAL(output)[output_ndx] = hSlab_buffer[ndx];
			}
			
		} else {
			//
			// no frames (i.e. a 3-d volume)
			if ( R_DEBUG_mincIO ) Rprintf("Debug: About to read value in 3-D volume\n");
			result = miget_real_value(minc_volume, &hSlab_start[1], 3, hSlab_buffer);
			if ( result != MI_NOERROR ) {
				error("Error in miget_real_value.  File: %s.\n", CHAR(STRING_ELT(filenames, i)));
			}
			REAL(output)[i] = *hSlab_buffer;
			
		}

		// done with this volume, so close it
		miclose_volume(minc_volume);
	}


	// // allocate the hyper-slab output buffer to the heap  (no need when using R_alloc)
	// if ( R_DEBUG_mincIO ) Rprintf("Freeing output buffer on heap.\n");
	// free(hSlab_buffer);

	// remove protection and return matrix
	UNPROTECT(1);
	if ( R_DEBUG_mincIO ) Rprintf("read_voxel_from_files: returning ...\n");
	return(output);
}





/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * Purpose: 
 *          
 *          
 *          
 *
 *          
*/
SEXP read_hyperslab(SEXP filename,  SEXP start,  SEXP count, SEXP nDimensions) {

	mihandle_t		minc_volume;
	int				result;
	int 			ndx;
	int				no_dimensions;


	SEXP			hSlab_buffer;
	double			*hSlab_buffer_ptr;
	unsigned long	hSlab_start[MI2_MAX_VAR_DIMS];
	unsigned long	hSlab_count[MI2_MAX_VAR_DIMS];
	int				hSlab_buffer_size;
	static char *dimorder3d[] = { "zspace","yspace","xspace" };
	static char *dimorder4d[] = { "time", "zspace","yspace","xspace" };



	// start ...
	if ( R_DEBUG_mincIO ) Rprintf("read_hyperslab: start ...\n");



	/* init number of dimensions in input volume
	  ... yes, I could get this myself, but it's best set in the calling R code */
	no_dimensions = INTEGER(nDimensions)[0];


	// init the hSlab start and count vectors
	hSlab_buffer_size = 1;
	for ( ndx=0; ndx < no_dimensions; ++ndx) {
		hSlab_start[ndx] = (unsigned long) INTEGER(start)[ndx];
		hSlab_count[ndx] = (unsigned long) INTEGER(count)[ndx];
		hSlab_buffer_size = hSlab_buffer_size * hSlab_count[ndx];
	}

	// allocate receiving matrix (3-D volumes always only return 1 column)
	if ( R_DEBUG_mincIO ) Rprintf("allocVector: size is %d doubles\n", hSlab_buffer_size);
	PROTECT(hSlab_buffer=allocVector(REALSXP, hSlab_buffer_size));
	hSlab_buffer_ptr = REAL(hSlab_buffer);


	// open the input volume
	result = miopen_volume(CHAR(STRING_ELT(filename,0)), MI2_OPEN_READ, &minc_volume);
	//
	if (result != MI_NOERROR) {
		error("Error opening input file: %s.\n", CHAR(STRING_ELT(filename, 0)));
	}

	/* set the apparent order to something conventional */
	if ( R_DEBUG_mincIO ) Rprintf("read_hyperslab: Setting the apparent order for %d dimensions ... ", no_dimensions);
	if ( no_dimensions == 3 ) {
		result = miset_apparent_dimension_order_by_name(minc_volume, 3, dimorder3d);
	} else if ( no_dimensions == 4 ) {
		result = miset_apparent_dimension_order_by_name(minc_volume, 4, dimorder4d);
	} else {
		error("Error file %s has %d dimensions and we can only deal with 3 or 4.\n", CHAR(STRING_ELT(filename,0)), no_dimensions);
	}
	if ( result != MI_NOERROR ) { 
		error("Error returned from miset_apparent_dimension_order_by_name while setting apparent order for %d dimensions.\n", no_dimensions); 
	}
	if ( R_DEBUG_mincIO ) Rprintf("Done.\n");

	// read the hyperslab
	result = miget_real_value_hyperslab(minc_volume, MI_TYPE_DOUBLE, 
														hSlab_start,
														hSlab_count,
														hSlab_buffer_ptr);
	if ( result != MI_NOERROR ) {
		error("Error in miget_real_value_hyperslab: %s.\n", CHAR(STRING_ELT(filename, 0)));
	}
	


	// remove protection and return matrix
	UNPROTECT(1);
	if ( R_DEBUG_mincIO ) Rprintf("read_hyperslab: returning ...\n");
	return(hSlab_buffer);
}






