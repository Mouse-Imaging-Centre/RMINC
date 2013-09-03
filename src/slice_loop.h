#ifndef __SLICE_LOOP_H__
#define __SLICE_LOOP_H__

#include <minc2.h>
#include <stdio.h>

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Utils.h>

extern mihandle_t* open_minc_files(SEXP filenames, 
				      misize_t *sizes);

// get_mask is broken - see per_voxel_anova for how to work around it
// for the moment
extern void get_mask(SEXP filename, 
		     mihandle_t hmask, 
		     double *mask_buffer,
		     misize_t *sizes);

extern double** create_slice_buffer(SEXP filenames,
				    misize_t *sizes);

extern void get_indices(int v0, int v1, int v2, 
			misize_t *sizes,
			misize_t *output_index, 
			misize_t *buffer_index);

extern void fill_slice_buffer(SEXP filenames, 
			      misize_t *sizes,
			      mihandle_t *hvol,
			      double **buffer,
			      int slice_number);

extern void get_mask_slice(mihandle_t hmask,
			   misize_t *sizes,
			   double *mask_buffer,
			   int slice_number);

extern void free_slice_buffer(SEXP filenames,
			      misize_t *sizes,
			      double **buffer);

extern void free_minc_files(SEXP filenames, 
			    mihandle_t *hvol);


#endif 
