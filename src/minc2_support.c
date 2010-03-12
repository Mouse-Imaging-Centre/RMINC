#include "minc2_io.h"


/*
   Given a minc filename, return a list containing:
   (1) the dimension names
   (2) the dimension sizes
   (3) and much, much more
   */
SEXP get_volume_info(SEXP filename) {

    mihandle_t          minc_volume;
	midimhandle_t		*dimensions;
	miclass_t			volume_class;
	mitype_t			volume_type;
	
	int					result, i;
	int					n_dimensions;
	int					n_protects, list_index;
	int					n_frames;
	
// variables to hold dim-related info
	unsigned int		dim_sizes[MI2_MAX_VAR_DIMS];
	double				dim_starts[MI2_MAX_VAR_DIMS];
	double				dim_steps[MI2_MAX_VAR_DIMS];
	double				time_offsets[MAX_FRAMES];
	double				time_widths[MAX_FRAMES];
	char 				*dim_name;
	char 				*dim_units;
	char 				*space_type;
	Rboolean			time_dim_exists;
	static char *dimorder3d[] = { "zspace","yspace","xspace" };
	static char *dimorder4d[] = { "time", "zspace","yspace","xspace" };
	

	/* declare R datatypes  */
	SEXP rtnList, listNames;
	SEXP xDimSizes, xDimNames, xDimUnits, xDimStarts, xDimSteps, xTimeWidths, xTimeOffsets;



	// start ...
	if ( R_DEBUG_mincIO ) Rprintf("get_volume_info: start ...\n");


	/* do some initialization */
	for (i=0; i < MI2_MAX_VAR_DIMS; ++i){					// set dim info to zeros
		dim_sizes[i] = 0;
		dim_starts[i] = 0;
		dim_steps[i] = 0;
	}

	// frame-related init
	time_dim_exists = FALSE;
	for (i=0; i < MAX_FRAMES; ++i) {
		time_offsets[i]=999.9;
		time_widths[i]=999.9;
	}
	n_frames = 0;

	n_protects = 0;								// counter of protected R variables



	/* init the return list (include list names) */
	PROTECT(rtnList=allocVector(VECSXP, R_RTN_LIST_LEN));
	PROTECT(listNames=allocVector(STRSXP, R_RTN_LIST_LEN));
	n_protects = n_protects +2;


	/* open the existing volume */
	result = miopen_volume(CHAR(STRING_ELT(filename,0)), MI2_OPEN_READ, &minc_volume);
	/* error on open? */
	if (result != MI_NOERROR) {
		error("Error opening input file: %s.\n", CHAR(STRING_ELT(filename,0)));
	}

	/* set the apparent order to something conventional */
	//	... first need to get the number of dimensions
	if ( miget_volume_dimension_count(minc_volume, MI_DIMCLASS_ANY, MI_DIMATTR_ALL, &n_dimensions) != MI_NOERROR ){
		error("Error returned from miget_volume_dimension_count.\n");
	}
	// ... now set the order
	if ( R_DEBUG_mincIO ) Rprintf("Setting the apparent order for %d dimensions ... ", n_dimensions);
	if ( n_dimensions == 3 ) {
		result = miset_apparent_dimension_order_by_name(minc_volume, 3, dimorder3d);
	} else if ( n_dimensions == 4 ) {
		result = miset_apparent_dimension_order_by_name(minc_volume, 4, dimorder4d);
	} else {
		error("Error file %s has %d dimensions and we can only deal with 3 or 4.\n", CHAR(STRING_ELT(filename,0)), n_dimensions);
	}
	if ( result != MI_NOERROR ) { 
		error("Error returned from miset_apparent_dimension_order_by_name while setting apparent order for %d dimensions.\n", n_dimensions); 
	}
	if ( R_DEBUG_mincIO ) Rprintf("Done.\n");

	/* get the volume data class (the intended "real" values) */ 
	if ( miget_data_class(minc_volume, &volume_class) != MI_NOERROR ){
		error("Error returned from miget_data_class.\n");
	}
	/* append to return list ... */
	list_index = 0;
	SET_VECTOR_ELT(rtnList, list_index, ScalarInteger(volume_class));
	SET_STRING_ELT(listNames, list_index, mkChar("volumeDataClass"));


	/* print the volume data type (as it is actually stored in the volume) */
	if ( miget_data_type(minc_volume, &volume_type) != MI_NOERROR ){
		error("Error returned from miget_data_type.\n");
	}
	/* append to return list ... */
	list_index++;
	SET_VECTOR_ELT(rtnList, list_index, ScalarInteger(volume_type));
	SET_STRING_ELT(listNames, list_index, mkChar("volumeDataType"));


	/* retrieve the volume space type (talairach, native, etc) */
	result = miget_space_name(minc_volume, &space_type);
	if ( result == MI_NOERROR ) { error("Error returned from miget_space_name.\n"); }
	/* append to return list ... */
	list_index++;
	SET_VECTOR_ELT(rtnList, list_index, mkString(space_type));
	SET_STRING_ELT(listNames, list_index, mkChar("spaceType"));


	/* retrieve the total number of dimensions in this volume */
	if ( miget_volume_dimension_count(minc_volume, MI_DIMCLASS_ANY, MI_DIMATTR_ALL, &n_dimensions) != MI_NOERROR ){
		error("Error returned from miget_volume_dimension_count.\n");
	}
	/* append to return list ... */
	list_index++;
	SET_VECTOR_ELT(rtnList, list_index, ScalarInteger(n_dimensions));
	SET_STRING_ELT(listNames, list_index, mkChar("nDimensions"));


	/* load up dimension-related information */
	//
	/* first allocate the R variables */
	PROTECT( xDimSizes=allocVector(INTSXP,n_dimensions) );
	PROTECT( xDimNames=allocVector(STRSXP,n_dimensions) );
	PROTECT( xDimUnits=allocVector(STRSXP,n_dimensions) );
	PROTECT( xDimStarts=allocVector(REALSXP,n_dimensions) );
	PROTECT( xDimSteps=allocVector(REALSXP,n_dimensions) );
	n_protects = n_protects +5;

	/* next, load up the midimension struct for all dimensions*/
	dimensions = (midimhandle_t *) malloc( sizeof( midimhandle_t ) * n_dimensions );
	result = miget_volume_dimensions(minc_volume, MI_DIMCLASS_ANY, MI_DIMATTR_ALL, MI_DIMORDER_APPARENT, n_dimensions, dimensions);
	// need to check against MI_ERROR, as "result" will contain nDimensions if OK
	if ( result == MI_ERROR ) { error("Error code(%d) returned from miget_volume_dimensions.\n", result); }


	/* get the dimension sizes for all dimensions */
	result = miget_dimension_sizes(dimensions, n_dimensions, dim_sizes);
	if ( result != MI_NOERROR ) { error("Error returned from miget_dimension_sizes.\n"); }
	/* add to R vector ... */
	for (i=0; i<n_dimensions; ++i){
		INTEGER(xDimSizes)[i] = dim_sizes[i];
	}
	list_index++;
	SET_VECTOR_ELT(rtnList, list_index, xDimSizes);
	SET_STRING_ELT(listNames, list_index, mkChar("dimSizes"));


	/* get the dimension START values for all dimensions */
	result = miget_dimension_starts(dimensions, MI_ORDER_FILE, n_dimensions, dim_starts);
	if ( result == MI_ERROR ) { error("Error returned from miget_dimension_starts.\n"); }
	/* add to R vector ... */
	for (i=0; i<n_dimensions; ++i){
		REAL(xDimStarts)[i] = dim_starts[i];
	}
	list_index++;
	SET_VECTOR_ELT(rtnList, list_index, xDimStarts);
	SET_STRING_ELT(listNames, list_index, mkChar("dimStarts"));


	/* get the dimension STEP values for all dimensions */
	result = miget_dimension_separations(dimensions, MI_ORDER_FILE, n_dimensions, dim_steps);
	if ( result == MI_ERROR ) { error("Error returned from miget_dimension_separations.\n"); }
	/* add to R vector ... */
	for (i=0; i<n_dimensions; ++i){
		REAL(xDimSteps)[i] = dim_steps[i];
	}
	list_index++;
	SET_VECTOR_ELT(rtnList, list_index, xDimSteps);
	SET_STRING_ELT(listNames, list_index, mkChar("dimSteps"));


	/* Loop over the dimensions to grab the remaining info ... */
	for( i=0; i < n_dimensions; ++i ){
	//
	/* get (and print) the dimension names for all dimensions*
	... remember that since miget_dimension_name calls strdup which, in turn,
	... calls malloc to get memory for the new string -- we need to call "mifree" on
	... our pointer to release that memory.  */
		result = miget_dimension_name(dimensions[i], &dim_name);
		
		// do we have a time dimension?
		if ( !strcmp(dim_name, "time") ) { 
			time_dim_exists = TRUE;
			n_frames = ( time_dim_exists ) ? dim_sizes[0] : 0;
		}
		
		// store the dimension name and units
		SET_STRING_ELT(xDimNames, i, mkChar(dim_name));
		mifree_name(dim_name);
		
		result = miget_dimension_units(dimensions[i], &dim_units);
		SET_STRING_ELT(xDimUnits, i, mkChar(dim_units));
		mifree_name(dim_units);
		
	}
	/* add number of frames to return list */
	list_index++;
	SET_VECTOR_ELT(rtnList, list_index, ScalarInteger(n_frames));
	SET_STRING_ELT(listNames, list_index, mkChar("nFrames"));
	
	// add dim names to return list
	list_index++;
	SET_VECTOR_ELT(rtnList, list_index, xDimNames);
	SET_STRING_ELT(listNames, list_index, mkChar("dimNames"));
	// add dim units
	list_index++;
	SET_VECTOR_ELT(rtnList, list_index, xDimUnits);
	SET_STRING_ELT(listNames, list_index, mkChar("dimUnits"));


	/* get the dimension OFFSETS values for the TIME dimension */
	if ( time_dim_exists ) {

		PROTECT( xTimeOffsets=allocVector(REALSXP,n_frames) );
		n_protects++;
		result = miget_dimension_offsets(dimensions[0], n_frames, 0, time_offsets);
		if ( result == MI_ERROR ) { error("Error returned from miget_dimension_offsets.\n"); }
		/* add to R vector ... */
		for (i=0; i < n_frames; ++i) {
			REAL(xTimeOffsets)[i] = time_offsets[i];
//			if (R_DEBUG_mincIO) Rprintf("Time offset[%d] =  %g\n", i, time_offsets[i]);
		}
		list_index++;
		SET_VECTOR_ELT(rtnList, list_index, xTimeOffsets);
		SET_STRING_ELT(listNames, list_index, mkChar("timeOffsets"));

		/* get the dimension WIDTH values for the TIME dimension */
		PROTECT( xTimeWidths=allocVector(REALSXP,n_frames) );
		n_protects++;
	
		result = miget_dimension_widths(dimensions[0], MI_ORDER_FILE, n_frames, 0, time_widths);
		if ( result == MI_ERROR ) { error("Error returned from miget_dimension_widths.\n"); }
		/* add to R vector ... */
		for (i=0; i<n_frames; ++i) {
			REAL(xTimeWidths)[i] = time_widths[i];
//			if (R_DEBUG_mincIO) Rprintf("Time width[%d] =  %g\n", i, time_widths[i]);
		}
		list_index++;
		SET_VECTOR_ELT(rtnList, list_index, xTimeWidths);
		SET_STRING_ELT(listNames, list_index, mkChar("timeWidths"));
	}



	// free heap memory
	free(dimensions);


	/* close volume */
	miclose_volume(minc_volume);


	/* attach the list component names to the list */
	setAttrib(rtnList, R_NamesSymbol, listNames);


	/* remove gc collection protection */
	UNPROTECT(n_protects);

   /* return */
	if ( R_DEBUG_mincIO ) Rprintf("get_volume_info: returning ...\n");
   return(rtnList);
}




/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 *  Purpose:
 *
 *    Convert voxel coordinates to a vector of world coords
 *
 * NOTE: 
 *    (1) It is assumed that we are receiving 0-relative voxel coordinates
 *    (2) input voxel coordinates are in volume order
 *    (3) returned world coordinates are in xyz order
 *
 * * * * * * * * * * * * * */
SEXP convert_voxel_to_world_mincIO(SEXP filename, SEXP voxCoords) {

	int			result, ndx;
	int			n_dimensions, offset;
	mihandle_t	minc_volume;
	double		voxel_coords[MI2_MAX_VAR_DIMS];
	double		world_coords[MI2_3D];
	static char *dimorder3d[] = { "zspace","yspace","xspace" };
	static char *dimorder4d[] = { "time", "zspace","yspace","xspace" };
	
	SEXP		output;

	// start ...
	if ( R_DEBUG_mincIO ) Rprintf("convert_voxel_to_world: start ...\n");



	/* open the volume */
	result = miopen_volume(CHAR(STRING_ELT(filename,0)), MI2_OPEN_READ, &minc_volume);
	//
	if (result != MI_NOERROR) {
		error("Error opening input file: %s.\n", CHAR(STRING_ELT(filename, 0)));
	}

	/* set the apparent order to something conventional */
	//	... first need to get the number of dimensions
	if ( miget_volume_dimension_count(minc_volume, MI_DIMCLASS_ANY, MI_DIMATTR_ALL, &n_dimensions) != MI_NOERROR ){
		error("Error returned from miget_volume_dimension_count.\n");
	}
	/* ... now set the order */
	if ( R_DEBUG_mincIO ) Rprintf("convert_world_to_voxel: Setting the apparent order for %d dimensions ... ", n_dimensions);
	if ( n_dimensions == 3 ) {
		result = miset_apparent_dimension_order_by_name(minc_volume, 3, dimorder3d);
	} else if ( n_dimensions == 4 ) {
		result = miset_apparent_dimension_order_by_name(minc_volume, 4, dimorder4d);
	} else {
		error("Error file %s has %d dimensions and we can only deal with 3 or 4.\n", CHAR(STRING_ELT(filename,0)), n_dimensions);
	}
	if ( result != MI_NOERROR ) { 
		error("Error returned from miset_apparent_dimension_order_by_name while setting apparent order for %d dimensions.\n", n_dimensions); 
	}
	if ( R_DEBUG_mincIO ) Rprintf("Done.\n");


	// init the voxel coordinates buffer
	for ( ndx=0; ndx < MI2_MAX_VAR_DIMS; ++ndx ) {
		voxel_coords[ndx] = 0;
	}

	offset = ( n_dimensions == 3 ) ? 0 : 1;
	for ( ndx=0; ndx < MI2_3D; ++ndx ) {
		voxel_coords[ndx +offset] = REAL(voxCoords)[ndx];
	}


	// do the call
	miconvert_voxel_to_world(minc_volume, voxel_coords, world_coords);


	// allocate output vector and then copy the results into it
	PROTECT(output=allocVector(REALSXP, MI2_3D));
	for ( ndx=0; ndx < MI2_3D; ++ndx ) {
		REAL(output)[ndx] = world_coords[ndx];
	}


	// remove protection and return vector
	UNPROTECT(1);
	if ( R_DEBUG_mincIO ) Rprintf("convert_voxel_to_world: returning ...\n");
	return(output);
}




/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 *  Purpose:
 *
 *    Convert world coordinates to a vector of voxel coords
 *
 * NOTE: 
 *    (1) 0-relative voxel coordinates are being returned
 *    (2) output voxel coordinates are in volume order, which we are forcing to be (t)zyx
 *    (3) input world coordinates are in xyz order
 *    (4) if the volume has 4D, then 4 voxel coords are returned, else 3 of 'em
 *
 * * * * * * * * * * * * * */
SEXP convert_world_to_voxel_mincIO(SEXP filename, SEXP worldCoords) {

	int			result, ndx;
	int			n_dimensions;
	mihandle_t	minc_volume;
	double		voxel_coords[MI2_MAX_VAR_DIMS];
	double		world_coords[MI2_3D];
	static char *dimorder3d[] = { "zspace","yspace","xspace" };
	static char *dimorder4d[] = { "time", "zspace","yspace","xspace" };
	
	SEXP		output;



	// start ...
	if ( R_DEBUG_mincIO ) Rprintf("convert_world_to_voxel: start ...\n");

	// init the world coordinates
	for ( ndx=0; ndx < MI2_3D; ++ndx ) {
		world_coords[ndx] = REAL(worldCoords)[ndx];
	}

	/* open the volume */
	result = miopen_volume(CHAR(STRING_ELT(filename,0)), MI2_OPEN_READ, &minc_volume);
	//
	if (result != MI_NOERROR) {
		error("Error opening input file: %s.\n", CHAR(STRING_ELT(filename, 0)));
	}


	/* set the apparent order to something conventional */
	//	... first need to get the number of dimensions
	if ( miget_volume_dimension_count(minc_volume, MI_DIMCLASS_ANY, MI_DIMATTR_ALL, &n_dimensions) != MI_NOERROR ){
		error("Error returned from miget_volume_dimension_count.\n");
	}
	/* ... now set the order */
	if ( R_DEBUG_mincIO ) Rprintf("convert_world_to_voxel: Setting the apparent order for %d dimensions ... ", n_dimensions);
	if ( n_dimensions == 3 ) {
		result = miset_apparent_dimension_order_by_name(minc_volume, 3, dimorder3d);
	} else if ( n_dimensions == 4 ) {
		result = miset_apparent_dimension_order_by_name(minc_volume, 4, dimorder4d);
	} else {
		error("Error file %s has %d dimensions and we can only deal with 3 or 4.\n", CHAR(STRING_ELT(filename,0)), n_dimensions);
	}
	
	if ( result != MI_NOERROR ) { 
		error("Error returned from miset_apparent_dimension_order_by_name while setting apparent order for %d dimensions.\n", n_dimensions); 
	}
	if ( R_DEBUG_mincIO ) Rprintf("Done.\n");


	// do the call
	miconvert_world_to_voxel(minc_volume,  world_coords, voxel_coords);

	// allocate output vector and then copy the results into it
	PROTECT(output=allocVector(REALSXP, n_dimensions));
	for ( ndx=0; ndx < n_dimensions; ++ndx ) {
		REAL(output)[ndx] = voxel_coords[ndx];
	}

	// remove protection and return vector
	UNPROTECT(1);
	if ( R_DEBUG_mincIO ) Rprintf("convert_world_to_voxel: returning ...\n");
	return(output);
}












