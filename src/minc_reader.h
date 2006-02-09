#ifndef __MINC_READER_H__
#define __MINC_READER_H__

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Utils.h>

/* function definitions */
extern void get_volume_sizes(char **filename, unsigned int *sizes);
extern void get_voxel_from_files(char **filenames, int *num_files,
				 int *v1, int *v2, int *v3, double *voxel);
extern void get_hyperslab(char **filename, int *start, int *count, 
			  double *slab);




#endif
