#include <Rcpp.h>
#include "minc2.h"
#include "minc_cpp.h"
#include <sstream>
using namespace Rcpp;
using namespace std;

typedef pair<RObject, int> indexed_robj;
bool comparator(const indexed_robj& l, const indexed_robj& r){
  return l.second < r.second;
}

// [[Rcpp::export]]
List rcpp_minc_apply(CharacterVector filenames,
                     bool use_mask,
                     CharacterVector mask,
                     double mask_lower_val,
                     double mask_upper_val,
                     RObject value_for_mask,
                     bool filter_masked,
                     NumericVector slab_sizes,
                     bool return_indices,
                     Function fun, 
                     List args) {

  vector<mihandle_t> volumes = open_minc2_volumes(filenames);
  mihandle_t mask_handle;
  vector<mihandle_t>::iterator volume_iterator;
  
  if(!check_same_dimensions(volumes)){
    throw range_error("At least one file is a different size");
  }
  
  if(use_mask){
    mask_handle = open_minc2_volumes(mask)[0];
    
    vector<mihandle_t> mask_and_vol;
    mask_and_vol.push_back(mask_handle);
    mask_and_vol.push_back(volumes[0]);
    
    if(!check_same_dimensions(mask_and_vol)){
      throw range_error("The mask and files differ in size");
    }  
  }
  
  midimhandle_t dimensions[3];
  misize_t sizes[3];  
  miget_volume_dimensions(volumes[0], MI_DIMCLASS_SPATIAL,
                          MI_DIMATTR_ALL, MI_DIMORDER_FILE,
                          3, dimensions);
  
  miget_dimension_sizes( dimensions, 3, sizes);
  
  misize_t hyperslab_dims[3];
  misize_t n_slabs[3];
  
  for(int i = 0; i < 3; ++i){
    hyperslab_dims[i] = (misize_t) slab_sizes[i];
    if(sizes[i] % hyperslab_dims[i] != 0){
      throw range_error("Volume not an even multiple of hyperslab size");
    }
    n_slabs[i] = sizes[i] / hyperslab_dims[i];
  }

  int nvols = volumes.size();
  NumericVector voxel_values(nvols);
  misize_t voxel_offsets[3];
  List output = List::create();
  double mask_buffer[hyperslab_dims[0]][hyperslab_dims[1]][hyperslab_dims[2]];
  double slab_buffer[nvols][hyperslab_dims[0]][hyperslab_dims[1]][hyperslab_dims[2]];
  int voxel_pos = 0;
  indexed_robj res;
  vector<indexed_robj> results;
  
  for(misize_t i = 0; i < n_slabs[0]; ++i){
    for(misize_t j= 0; j < n_slabs[1]; ++j){
      for(misize_t k = 0; k < n_slabs[2]; ++k ){
        
        //Set voxel coords
        voxel_offsets[0] = i * hyperslab_dims[0];
        voxel_offsets[1] = j * hyperslab_dims[1];
        voxel_offsets[2] = k * hyperslab_dims[2];
  
        
        //Check if a mask was supplied, then check if the current voxel is masked
        if(use_mask){
          cautious_get_hyperslab(mask_handle,    // read from handle
                                 MI_TYPE_DOUBLE, // double data
                                 voxel_offsets,  // starting from position
                                 hyperslab_dims, // how many voxels
                                 &mask_buffer,   // into
                                 "Error Reading Mask");  
        }
        
        for(int vol = 0; vol < nvols; ++vol){
          stringstream error_message;
          error_message.str("Error Reading Volume ");
          error_message << (vol + 1);
          
          cautious_get_hyperslab(volumes[vol],
                                 MI_TYPE_DOUBLE,
                                 voxel_offsets,
                                 hyperslab_dims,
                                 &slab_buffer[vol][0][0][0],
                                 error_message.str());
        }
        
        for(misize_t x = 0; x < hyperslab_dims[0]; ++x){
          for(misize_t y = 0; y < hyperslab_dims[1]; ++y){
            for(misize_t z = 0; z < hyperslab_dims[2]; ++z){
              for(int xvol = 0; xvol < nvols; ++xvol){
                voxel_values[xvol] = slab_buffer[xvol][x][y][z];
              }
              
              voxel_pos = (int)(sizes[1] * sizes[2] * (voxel_offsets[0] + x) + 
                sizes[2] * (voxel_offsets[1] + y) + 
                (voxel_offsets[2] + z));

              if((!use_mask) || 
                 (mask_buffer[x][y][z] > (mask_lower_val - .5) &&
                 mask_buffer[x][y][z] < (mask_upper_val + .5))){
                
                res = make_pair(fun(voxel_values, args), voxel_pos);
                results.push_back(res);
                
              } else if(!filter_masked){
                res = make_pair(value_for_mask, voxel_pos);
                results.push_back(res);
              }
            }
          }
        }
      }
    }
  }
  
  sort(results.begin(), results.end(), comparator);
  
  vector<indexed_robj>::iterator index_peeler;
  for(index_peeler = results.begin(); index_peeler != results.end(); ++index_peeler){
    indexed_robj current_res = *index_peeler;
    
    if(return_indices){
      List indexed_res = List::create(current_res.first, current_res.second);
      output.push_back(indexed_res);
    } else {
      output.push_back(current_res.first);
    }
  }
  
  for(volume_iterator = volumes.begin(); volume_iterator != volumes.end(); ++volume_iterator){
    miclose_volume(*volume_iterator);
  }
  
  if(use_mask){
    miclose_volume(mask_handle);
  }
  
  return(output);
}
