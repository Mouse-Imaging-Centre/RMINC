#ifndef __MINC_APPLY_H__
#define __MINC_APPLY_H__

#include <Rcpp.h>
#include "minc2.h"
using namespace Rcpp;
using namespace std;


extern mihandle_t open_minc2_volume(CharacterVector filename);
extern vector<mihandle_t> open_minc2_volumes(CharacterVector filenames);
extern bool check_same_dimensions(vector<mihandle_t> volumes);
extern List rcpp_minc_apply(CharacterVector filenames,
                            bool use_mask,
                            CharacterVector mask,
                            double mask_lower_val,
                            double mask_upper_val,
                            RObject value_for_mask,
                            bool filter_masked,
                            Function fun, 
                            List args);

extern void cautious_get_hyperslab(mihandle_t volume,
                                   mitype_t buffer_data_type,
                                   misize_t *voxel_offsets,
                                   misize_t *sizes,
                                   void *buffer,
                                   String error_message);

extern void cautious_open_volume(char *filename, 
                                 int mode, 
                                 mihandle_t *volume, 
                                 String error_message)

#endif