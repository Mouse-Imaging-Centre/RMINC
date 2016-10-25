#include <Rcpp.h>
#include "minc2.h"
#include <stack>
#include "minc_cpp.h"
using namespace Rcpp;
using namespace std;

template <typename T>
string NumberToString ( T Number )
{
  ostringstream ss;
  ss.clear();
  ss << Number;
  return ss.str();
}

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
  vector<misize_t> volume_dimensions(3, 0);
  
  int success = miget_volume_dimensions(volume, MI_DIMCLASS_SPATIAL,
                                        MI_DIMATTR_ALL, MI_DIMORDER_FILE,
                                        3, dimensions);
  
  if(success == MI_ERROR){
    stop("Couldn't read volume dimensions");
  }
  
  success = miget_dimension_sizes(dimensions, 3, sizes);
  
  if(success != MI_NOERROR){
    stop("Couldn't read dimension sizes");
  }
  
  for(int i = 0; i < 3; ++i){
    volume_dimensions[i] = sizes[i];
  }
  
  return(volume_dimensions);
}

CharacterVector path_to_filename(CharacterVector filenames){
  Environment base("package:base");
  Function strsplit = base["strsplit"];
  Function paste = base["paste"];
  
  List split_res(filenames.size());
  CharacterVector lasts(filenames.size());
  NumericVector lengths(filenames.size());
  for(int i = 0; i < filenames.size(); ++i){
    CharacterVector res = as<CharacterVector>(as<List>(strsplit(as<std::string>(filenames[i]), "/"))[0]);

    for(int r = 0; r < res.size(); ++r) //Remove empty substrings
      if(res[r] == "") res.erase(res.begin() + r, res.begin() + r + 1);
    
    split_res[i] = res;
    lengths[i] = res.size();
    lasts[i] = res[lengths[i] - 1];
  }
  
  // Check if there are any duplicates in last
  if(unique(lasts).size() == lasts.size())
    return(lasts);
  
  int longest_possible_sub = (int) min(lengths) - 1;
  bool all_equal = true;
  for(int i = 0; i < longest_possible_sub; ++i){
    for(int f = 1; f < filenames.size(); ++f){
      if(as<CharacterVector>(split_res[0])[0] != as<CharacterVector>(split_res[f])[0]){
        all_equal = false;
        break;
      }
    }
    
    // if they\'re all the same trim them
    if(all_equal){
      for(int f = 0; f < filenames.size(); ++f){
        CharacterVector res = split_res[f];
        CharacterVector::iterator it = res.begin();
        res.erase(it, it + 1);
        split_res[f] = res;
      }
    }
  }
  
  CharacterVector results(filenames.size());
  for(int f = 0; f < filenames.size(); ++f){
    results[f] = as<CharacterVector>(paste(split_res[f], Named("collapse", "_")))[0];
  }
  
  return(results); 
}