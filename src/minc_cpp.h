#ifndef __MINC_APPLY_H__
#define __MINC_APPLY_H__

#include <Rcpp.h>
#include <numeric>
#include <memory>
#include <functional>
#include <set>
#include <stack>
#include "minc2.h"
using namespace Rcpp;
using namespace std;


extern SEXP get_volume(std::string filename);
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
extern vector<double> get_step_sizes(mihandle_t volume);

class MincVolume {
private:
  string filename;
  shared_ptr<mivolume> handle;
  vector<double> steps;
  vector<misize_t> sizes;
  bool debug;

  void init(string fn, int mode, bool dbg){
    filename = fn;
    mihandle_t hndl;
    int res = miopen_volume(fn.c_str(), mode, &hndl);
    if(res != MI_NOERROR){
      stop("unable to open file: " + filename + "\n");
    }

    debug = dbg;
    if(debug)
      Rcout << "Assigning handle: " << hndl << "\n";
    
    handle = shared_ptr<mivolume>(hndl, [=](mihandle_t h){
        if(debug)
          Rcout << "Deleting handle: " << h << "\n";
        miclose_volume(h);
      });
    steps = get_step_sizes(handle.get());
    sizes = get_volume_dimensions(handle.get());
    
  }
  
public:
  MincVolume(){}
    
  MincVolume(string fn, int mode = MI2_OPEN_READ, bool dbg = false){
    init(fn, mode, dbg);
  }

  mihandle_t get_handle(){ return(handle.get()); }
  string get_filename(){ return(filename); }
  vector<double> get_steps(){ return(steps); }
  vector<misize_t> get_sizes(){ return(sizes); }
  int size(){
    int res = (int) std::accumulate(sizes.begin(), sizes.end(), 1, std::multiplies<misize_t>());
    return(res);
  }
  double voxel_volume(){
    double res = fabs(std::accumulate(steps.begin(), steps.end(), 1.0, std::multiplies<double>()));
    return(res);
  }
  shared_ptr<double> read_volume(mitype_t type);
  void read_volume_to_buffer(double* &buf, mitype_t type);
  shared_ptr<double> read_slab(vector<misize_t> start, vector<misize_t> count, mitype_t type);
  void read_slab_to_buffer(vector<misize_t> start, vector<misize_t> count
                           , mitype_t type, double* &buf);
};

extern map<int,double> atlasVol(MincVolume subject, double vox_vol);
extern map<int,double> atlasMeans(MincVolume subject, shared_ptr<int> atlas_buf);
extern map<int,double> atlasReduce(MincVolume subject, shared_ptr<int> atlas_buf
                                   , std::function<double(double,double)>);
extern List mergeCountMaps(vector<map<int,double> > maps);
extern List atlas_get_all(CharacterVector filenames, CharacterVector atlas
                          , std::string method);

#endif
