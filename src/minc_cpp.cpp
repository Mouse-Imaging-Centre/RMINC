#include <Rcpp.h>
#include "minc2.h"
#include "minc_cpp.h"
using namespace Rcpp;
using namespace std;

void cautious_get_hyperslab(mihandle_t volume,
                            mitype_t buffer_data_type,
                            misize_t *voxel_offsets,
                            misize_t *sizes,
                            void *buffer,
                            string error_message){
  int res = miget_real_value_hyperslab(volume, buffer_data_type, voxel_offsets, sizes, buffer);
  if(res != MI_NOERROR){
    stop(error_message);
  }
}

void cautious_open_volume(char *filename, int mode, mihandle_t *volume, string error_message){
  int res = miopen_volume(filename, mode, volume);
  if(res != MI_NOERROR){
    stop(error_message);
  }
}

mihandle_t open_minc2_volume(CharacterVector filename){
  mihandle_t current_handle;
  int read_result;
  
  cautious_open_volume(filename[0],
                       MI2_OPEN_READ, 
                       &current_handle,
                       "Trouble reading file: " + filename[0]); 
  
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

vector<misize_t> get_volume_dimensions(mihandle_t volume){
  midimhandle_t dimensions[3];
  misize_t sizes[3];
  vector<misize_t> volume_dimensions;
  
  int success = miget_volume_dimensions(volume, MI_DIMCLASS_SPATIAL,
                                        MI_DIMATTR_ALL, MI_DIMORDER_FILE,
                                        3, dimensions);
  
  if(success != MI_NOERROR){
    stop("Couldn't read volume dimensions");
  }
  
  success = miget_dimension_sizes(dimensions, 3, sizes);
  
  if(success != MI_NOERROR){
    stop("Couldn't read dimension sizes");
  }
  
  for(int i = 0; i < 3; ++i){
    volume_dimensions.push_back(sizes[i]);
  }
  
  return(volume_dimensions);
}
