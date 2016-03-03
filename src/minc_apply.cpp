#include <Rcpp.h>
#include "minc2.h"
using namespace Rcpp;
using namespace std;


mihandle_t open_minc2_volume(CharacterVector filename){
  mihandle_t current_handle;
  int read_result;
  
  read_result = miopen_volume(filename[0],
                              MI2_OPEN_READ, &current_handle); 
  
  if(read_result != MI_NOERROR){
    throw range_error("Trouble reading file: " + filename[0]);
  }
  
  return(current_handle);
}

vector<mihandle_t> open_minc2_volumes(CharacterVector filenames){
  
  vector<mihandle_t> volumes;
  mihandle_t current_handle;
  CharacterVector::iterator file_iterator;
  vector<mihandle_t>::iterator volume_iterator;
  
  for(file_iterator = filenames.begin();
      file_iterator != filenames.end();
      ++file_iterator){
    try {
      current_handle = open_minc2_volume(wrap(*file_iterator));
    } catch(...){
      for(volume_iterator = volumes.begin(); volume_iterator != volumes.end(); ++volume_iterator){
        miclose_volume(*volume_iterator);
      }
      throw;
    }
    
    volumes.push_back(current_handle);
  }
  
  return(volumes);
}



bool check_same_dimensions(vector<mihandle_t> volumes){
  vector<mihandle_t>::iterator volume_iterator;
  
  midimhandle_t first_dims[3];
  misize_t first_sizes[3];  
  midimhandle_t dimensions[3];
  misize_t sizes[3];
  
  miget_volume_dimensions(volumes[0], MI_DIMCLASS_SPATIAL,
                          MI_DIMATTR_ALL, MI_DIMORDER_FILE,
                          3, first_dims);
  
  miget_dimension_sizes( first_dims, 3, first_sizes);
  
  bool all_same_size = true;
  for(volume_iterator = volumes.begin() + 1; 
      volume_iterator != volumes.end() && all_same_size; 
      ++volume_iterator){
    
    miget_volume_dimensions(*volume_iterator, MI_DIMCLASS_SPATIAL,
                            MI_DIMATTR_ALL, MI_DIMORDER_FILE,
                            3, dimensions);
    
    miget_dimension_sizes(dimensions, 3, sizes);
    
    all_same_size = 
      all_same_size && 
      sizes[0] == first_sizes[0] &&
      sizes[1] == first_sizes[1] &&
      sizes[2] == first_sizes[2];
  }
  
  return(all_same_size);
}

// [[Rcpp::export]]
List rcpp_minc_apply(CharacterVector filenames,
                     bool use_mask,
                     CharacterVector mask,
                     double mask_lower_val,
                     double mask_upper_val,
                     RObject value_for_mask,
                     bool filter_masked,
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

  NumericVector voxel_values;
  double voxel_from_file;
  unsigned int nvols = volumes.size();
  misize_t voxel_coords[3];
  List output = List::create();
  bool is_masked = false;
  double mask_value;
  
  for(misize_t i = 0; i < sizes[0]; ++i){
    for(misize_t j= 0; j < sizes[1]; ++j){
      for(misize_t k = 0; k < sizes[2]; ++k ){
        
        //For each voxel ijk, re-initialize voxel_values
        voxel_values = NumericVector(nvols);
        
        //Set voxel coords
        voxel_coords[0] = i;
        voxel_coords[1] = j;
        voxel_coords[2] = k;
        
        is_masked = false;
        
        //Check if a mask was supplied, then check if the current voxel is masked
        if(use_mask){
          miget_real_value(mask_handle, voxel_coords, 3, &mask_value);
          
          if((mask_value < mask_lower_val - .5) ||
             (mask_value > mask_upper_val + .5)){
            is_masked = true;
          }
        }
        
        //If the voxel isn't masked, gather the data, run the function, otherwise push 
        if(!is_masked){
          //Read the values into voxel_values, 
          for(unsigned int vol = 0; vol < nvols; ++vol){
            miget_real_value(volumes[vol],
                             voxel_coords,
                             3,
                             &voxel_from_file);
            
            voxel_values[vol] = voxel_from_file;
          }
          //Run the user specified function
          output.push_back(fun(voxel_values, args));
        } else {
          //If the user wants to hang on to the masked values
          if(!filter_masked){
            output.push_back(value_for_mask); 
          }
        }
      }
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

/**
 * RMINC::minc_apply(c("/tmp/rminctestdata/absolute_jacobian_file_1.mnc", 
 * "/tmp/rminctestdata/absolute_jacobian_file_2.mnc", mean, list())
 */
