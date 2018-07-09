#ifndef __MINC_APPLY_H__
#define __MINC_APPLY_H__

#include <Rcpp.h>
#include <numeric>
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
                                   string error_message);

extern void cautious_open_volume(char *filename, 
                                 int mode, 
                                 mihandle_t *volume, 
                                 string error_message);

extern vector<misize_t> get_volume_dimensions(mihandle_t volume);

class MincVolume {
private:
  string filename;
  mihandle_t handle;
  vector<misize_t> sizes;

  void init(string fn, int mode){
    filename = fn;
    mihandle_t hndl;
    int res = miopen_volume(fn.c_str(), mode, &hndl);
    if(res != MI_NOERROR){
      stop("unable to open file: " + filename + "\n");
    }

    handle = hndl;
    sizes = get_volume_dimensions(handle);
  }
  
public:
  MincVolume(string fn, int mode){
    init(fn, mode);
  }

  MincVolume(string fn){
    init(fn, MI2_OPEN_READ);
  }
  
  ~MincVolume(){
    miclose_volume(handle);
  }
  
  mihandle_t get_handle(){ return(handle); }
  int size(){
    int res = (int) std::accumulate(sizes.begin(), sizes.end(), 1, std::multiplies<misize_t>());
    return(res);
  }
  double* read_volume(mitype_t type);
  void read_volume_to_buffer(double* &buf, mitype_t type);
  double* read_slab(vector<misize_t> start, vector<misize_t> count, mitype_t type);
  void read_slab_to_buffer(vector<misize_t> start, vector<misize_t> count
                           , mitype_t type, double* &buf);
};

#endif
