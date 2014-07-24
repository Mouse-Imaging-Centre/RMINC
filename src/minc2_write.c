#include "minc2_io.h"



/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * Purpose: 
 *          
 *          
 *          
 *
 *          
*/
SEXP write_volume(SEXP filename, SEXP nDimensions, 
								SEXP dimLengths, 
								SEXP dimStarts, 
								SEXP dimSteps, 
								SEXP volumeDataType, 
								SEXP volumeRange, 
								SEXP hSlab) {

	mihandle_t			minc_volume;
	midimhandle_t		dim[MI2_MAX_VAR_DIMS];
	mivolumeprops_t		volume_properties;
	
	int				result;
	int				ndx;
	int				no_dimensions;
	int				dim_lengths[MI2_MAX_VAR_DIMS];
	double			dim_starts[MI2_MAX_VAR_DIMS];
	double			dim_steps[MI2_MAX_VAR_DIMS];
	double			volume_range_min, volume_range_max;
	int				volume_data_type;
	
	// pointer to the volume data
	double			*hSlab_p;
	misize_t	hSlab_start[MI2_MAX_VAR_DIMS];
	misize_t	hSlab_count[MI2_MAX_VAR_DIMS];


	// start ...
	if ( R_DEBUG_mincIO ) Rprintf("write_volume: start ...\n");

	/* init number of dimensions and their respective sizes
	  ... yes, I could get this myself, but it's best set in the calling R code */
	no_dimensions = INTEGER(nDimensions)[0];
	volume_data_type = INTEGER(volumeDataType)[0];
	volume_range_min = REAL(volumeRange)[0];
	volume_range_max = REAL(volumeRange)[1];
	
	// init volume data pointer
	hSlab_p = REAL(hSlab);
	
	
	// init lenghts, steps, and starts
	for (ndx=0; ndx < no_dimensions; ++ndx) {
		dim_lengths[ndx]  = INTEGER(dimLengths)[ndx];
		hSlab_count[ndx]  = INTEGER(dimLengths)[ndx];
		dim_starts[ndx]  = REAL(dimStarts)[ndx];
		dim_steps[ndx]  = REAL(dimSteps)[ndx];
	}


	// set the properties for the new output volume
	// ... no compression, no chunking, no multi-resolution, ... nothin fancy
	result = minew_volume_props(&volume_properties);
	if (result != MI_NOERROR) {
		error("write_volume:minew_volume_props: Error setting output volume properties: %s.\n", CHAR(STRING_ELT(filename, 0)));
	}
	result = miset_props_compression_type(volume_properties, MI_COMPRESS_NONE);
	if (result != MI_NOERROR) {
		error("write_volume:miset_props_compression_type: Error setting output volume properties: %s.\n", CHAR(STRING_ELT(filename, 0)));
	}
	result = miset_props_multi_resolution(volume_properties, FALSE, 1);
	if (result != MI_NOERROR) {
		error("write_volume:miset_props_multi_resolution: Error setting output volume properties: %s.\n", CHAR(STRING_ELT(filename, 0)));
	}

	/*	create the 3 output dimensions in the order Z, Y, X, as the volume
		is stored in this order */
	// z-dim
	result = micreate_dimension("zspace", 
								MI_DIMCLASS_SPATIAL,
								MI_DIMATTR_REGULARLY_SAMPLED, 
								dim_lengths[0],
								&dim[0]);
	//
	if (result != MI_NOERROR) {
		error("write_volume: Error initializing the dimension struct for %s dimension.\n", "zspace");
	}
	
	// y-dim
	result = micreate_dimension("yspace", 
								MI_DIMCLASS_SPATIAL,
								MI_DIMATTR_REGULARLY_SAMPLED, 
								dim_lengths[1],
								&dim[1]);
	//
	if (result != MI_NOERROR) {
		error("write_volume: Error initializing the dimension struct for %s dimension.\n", "yspace");
	}
	
	// x-dim
	result = micreate_dimension("xspace", 
								MI_DIMCLASS_SPATIAL,
								MI_DIMATTR_REGULARLY_SAMPLED, 
								dim_lengths[2],
								&dim[2]);
	//
	if (result != MI_NOERROR) {
		error("write_volume: Error initializing the dimension struct for %s dimension.\n", "xspace");
	}


	// set the start values for each dimension
	result = miset_dimension_starts(dim, no_dimensions, dim_starts);
	if (result == MI_ERROR) {
		error("write_volume: Error setting dimension start values.\n");
	}

	// set the step values for each dimension
	result = miset_dimension_separations(dim, no_dimensions, dim_steps);
	if (result == MI_ERROR) {
		error("write_volume: Error setting dimension step values.\n");
	}


	// create the volume structure (no data yet, of course)
	result = micreate_volume(CHAR(STRING_ELT(filename, 0)),
								no_dimensions, 
								dim,
								volume_data_type,
								MI_CLASS_REAL,
								volume_properties,
								&minc_volume);
	//
	if (result != MI_NOERROR) {
		error("write_volume: Error creating output volume structure: %s.\n", CHAR(STRING_ELT(filename, 0)));
	}
		
	// create the path to the image data
	result = micreate_volume_image(minc_volume);
	//
	if (result != MI_NOERROR) {
		error("write_volume: Error writing data to volume %s.\n", CHAR(STRING_ELT(filename, 0)));
	}


	// set valid and real ranges 
	// ... 0xFFFF=65535=16-bits (unsigned)
	miset_volume_valid_range(minc_volume, 0x7FFF, 0);
	miset_volume_range(minc_volume, volume_range_max, volume_range_min);


	// write  hyperslab (entire volume)
	hSlab_start[0] = hSlab_start[1] = hSlab_start[2] = 0;
	if ( R_DEBUG_mincIO ) Rprintf("hSlab_count [0..2] = %d, %d, %d\n", 
		hSlab_count[0], hSlab_count[1], hSlab_count[2]);
			
	result = miset_real_value_hyperslab(minc_volume, MI_TYPE_DOUBLE, 
														hSlab_start,
														hSlab_count,
														hSlab_p);
	if ( result != MI_NOERROR ) {
		error("Error in miset_real_value_hyperslab: %s.\n", CHAR(STRING_ELT(filename, 0)));
	}




	// free resources
        //
        // in the current version of minc (libsrc2/volume.c), these volume properties
        // as well as the dimension handles are already freed by the miclose_volume function.
        //
	//mifree_volume_props(volume_properties);								// free the property list
	//for ( ndx=0; ndx<no_dimensions; ++ndx) {							// free the dimhandles
	//	mifree_dimension_handle(dim[ndx]);
	//}
	// close new volume
	result = miclose_volume(minc_volume);
	if (result != MI_NOERROR) {
		error("write_volume: Error closing newly created volume %s.\n", CHAR(STRING_ELT(filename, 0)));
	}

	// done, return NULL
	if ( R_DEBUG_mincIO ) Rprintf("write_volume: returning ...\n");
	return(R_NilValue);
}

