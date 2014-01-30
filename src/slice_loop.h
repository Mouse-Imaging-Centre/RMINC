#ifndef __SLICE_LOOP_H__
#define __SLICE_LOOP_H__

#include <minc2.h>
#include <stdio.h>

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Utils.h>

extern mihandle_t* open_minc_files(SEXP filenames, 
				      unsigned int *sizes);

// get_mask is broken - see per_voxel_anova for how to work around it
// for the moment
extern void get_mask(SEXP filename, 
		     mihandle_t hmask, 
		     double *mask_buffer,
		     unsigned int *sizes);

extern double** create_slice_buffer(SEXP filenames,
				    unsigned int *sizes);

extern void get_indices(int v0, int v1, int v2, 
			unsigned int *sizes,
			unsigned long *output_index, 
			unsigned long *buffer_index);

extern void fill_slice_buffer(SEXP filenames, 
			      unsigned int *sizes,
			      mihandle_t *hvol,
			      double **buffer,
			      int slice_number);

extern void get_mask_slice(mihandle_t hmask,
			   unsigned int *sizes,
			   double *mask_buffer,
			   int slice_number);

extern void free_slice_buffer(SEXP filenames,
			      unsigned int *sizes,
			      double **buffer);

extern void free_minc_files(SEXP filenames, 
			    mihandle_t *hvol);


#endif 
