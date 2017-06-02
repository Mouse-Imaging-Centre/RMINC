#include <Rcpp.h>
#include "minc2.h"
#include "minc_cpp.h"
#include <sstream>
#include <stdlib.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
List rcpp_minc_apply(CharacterVector filenames,
                     bool use_mask,
                     CharacterVector mask,
                     double mask_lower_val,
                     double mask_upper_val,
                     RObject value_for_mask,
                     bool filter_masked,
                     NumericVector slab_sizes,
                     Function fun, 
                     List args) {

  vector<mihandle_t> volumes = open_minc2_volumes(filenames);
  mihandle_t mask_handle;
  vector<mihandle_t>::iterator volume_iterator;
  
  if(!check_same_dimensions(volumes)){
    for(int i = 0; i < volumes.size(); ++i){
      miclose_volume(volumes[i]);
    }
    
    stop("At least one file is a different size");
  }
  
  if(use_mask){
    mask_handle = open_minc2_volumes(mask)[0];
    
    vector<mihandle_t> mask_and_vol;
    mask_and_vol.push_back(mask_handle);
    mask_and_vol.push_back(volumes[0]);
    
    if(!check_same_dimensions(mask_and_vol)){
      miclose_volume(mask_and_vol[0]);
      miclose_volume(mask_and_vol[1]);
      for(int i = 0; i < volumes.size(); ++i){
        miclose_volume(volumes[i]);
      }
      
      stop("The mask and files differ in size");
    }  
  }
  
  vector<misize_t> sizes = get_volume_dimensions(volumes[0]);
  misize_t hyperslab_dims[3];
  misize_t n_slabs[3];
  
  for(int i = 0; i < 3; ++i){
    hyperslab_dims[i] = (misize_t) slab_sizes[i];
    if(sizes[i] % hyperslab_dims[i] != 0){
      n_slabs[i] = 1;
    } else {
      n_slabs[i] = 0;
    }
    
    n_slabs[i] += (sizes[i] / hyperslab_dims[i]);
  }
  misize_t old_hyperslab_dims[3] = {hyperslab_dims[0], hyperslab_dims[1], hyperslab_dims[2]};

  int nvols = volumes.size();
  int nvoxels = (int) sizes[0] * sizes[1] * sizes[2];

  NumericVector voxel_values(nvols);
  misize_t voxel_offsets[3];

  //Setup buffers
  double *mask_buffer =
    (double *) calloc(hyperslab_dims[0] * hyperslab_dims[1] * hyperslab_dims[2], sizeof(double));
  
  double **slab_buffer =
    (double **) malloc(nvols * sizeof(double *));
  
  for(int i = 0; i < nvols; ++i){
    slab_buffer[i] = 
      (double *) malloc(hyperslab_dims[0] * hyperslab_dims[1] * hyperslab_dims[2] * sizeof(double));
  }

  // Setup looping constructs
  int voxel_pos = 0;
  int added_results = 0;
  List results(nvoxels);
  NumericVector result_inds(nvoxels);
  
  Rprintf("Number of Volumes: %d\n", nvols);
  Rprintf("Number of slabs: %d\n", n_slabs[0] * n_slabs[1] * n_slabs[2]);
  Rprintf("In slab: ");
 
 for(misize_t i = 0; i < n_slabs[0]; ++i){
   for(misize_t j= 0; j < n_slabs[1]; ++j){
     for(misize_t k = 0; k < n_slabs[2]; ++k ){
       
       Rprintf("%d ", (i * n_slabs[1] * n_slabs[2]) + 
         (j * n_slabs[2]) + 
         k);
       
       //Set voxel coords
       voxel_offsets[0] = i * hyperslab_dims[0];
       voxel_offsets[1] = j * hyperslab_dims[1];
       voxel_offsets[2] = k * hyperslab_dims[2];
       
       //Check for uneven final hyperslab
       for(int i = 0; i < 3; ++i){
         if((signed int)(sizes[i] - voxel_offsets[i] - hyperslab_dims[i]) < 0){
           hyperslab_dims[i] = sizes[i] - voxel_offsets[i];
         }
       }
       
       int hyperslab_vol = hyperslab_dims[0] * hyperslab_dims[1] * hyperslab_dims[2];
       
       //Check for user breaking the loop
       checkUserInterrupt();
       
       //Default to processing each slab
       bool process_slab = true;
       
       //Check if a mask was supplied, then check if the current voxel is masked
       if(use_mask){
         cautious_get_hyperslab(mask_handle,    // read from handle
                                MI_TYPE_DOUBLE, // double data
                                voxel_offsets,  // starting from position
                                hyperslab_dims, // how many voxels
                                mask_buffer,   // into
                                "Error Reading Mask");
         
         //If using a mask, only process as slab if atleast one voxel isn't masked
         process_slab = false;
         for(int mv = 0; mv < hyperslab_vol; ++mv){
           if(mask_buffer[mv] > (mask_lower_val - .5) &&
              mask_buffer[mv] < (mask_upper_val + .5)){
             process_slab = true;
           }
         }
       }
       
       
       //If the slab is to be processed, read in the subjects
       if(process_slab){
         for(int vol = 0; vol < nvols; ++vol){
           stringstream error_message;
           error_message << "Error Reading Volume " << (vol + 1) << "\n";
           
           cautious_get_hyperslab(volumes[vol],
                                  MI_TYPE_DOUBLE,
                                  voxel_offsets,
                                  hyperslab_dims,
                                  slab_buffer[vol],
                                             error_message.str());
           
           // Rprintf("Just read in a slab: ");
           // for(int i = 0; i < (hyperslab_dims[0] * hyperslab_dims[1] * hyperslab_dims[2]); ++i){
           //   Rprintf("%d, ", *(mask_buffer + i));
           // }
           // Rprintf("\n");
         }
       }
         
         for(misize_t x = 0; x < hyperslab_dims[0]; ++x){
           for(misize_t y = 0; y < hyperslab_dims[1]; ++y){
             for(misize_t z = 0; z < hyperslab_dims[2]; ++z){
               
               voxel_pos = (int)(sizes[1] * sizes[2] * (voxel_offsets[0] + x) +
                 sizes[2] * (voxel_offsets[1] + y) +
                 (voxel_offsets[2] + z));
               
               
               if(!process_slab && !filter_masked){
                 results[added_results] = value_for_mask;
                 result_inds[added_results] = voxel_pos;
                 ++added_results;
               } else {
                
                 //ugly indexing into 2nd-D of slab buffer
                 int slab_pos =  x * hyperslab_dims[1] * hyperslab_dims[2] +
                   y * hyperslab_dims[2] +
                   z;
                 
                 for(int xvol = 0; xvol < nvols; ++xvol){
                   voxel_values[xvol] = slab_buffer[xvol][slab_pos];
                 }
                 
                 double mask_val =
                   mask_buffer[slab_pos];
                 
                 if((!use_mask) ||
                    (mask_val > (mask_lower_val - .5) &&
                    mask_val < (mask_upper_val + .5))){
                   
                   //Rprintf("A few voxel values %d %d %d", voxel_values[0], voxel_values[1], voxel_values[2]);
                   results[added_results] = fun(voxel_values, args);
                   result_inds[added_results] = voxel_pos;
                   ++added_results;
                   
                 } else if(!filter_masked){
                   results[added_results] = value_for_mask;
                   result_inds[added_results] = voxel_pos;
                   ++added_results;
                 }
               }
             }
           }
         }
     
        
        //reset old hyperslab dims if they changed
        for(int i = 0; i < 3; ++i){
          hyperslab_dims[i] = old_hyperslab_dims[i];
        }
      }
    }
  }
  
  Rprintf("\n");

  for(int i = 0; i < nvols; ++i){
    free(slab_buffer[i]);
  }
  free(slab_buffer);
  free(mask_buffer);

  for(volume_iterator = volumes.begin(); volume_iterator != volumes.end(); ++volume_iterator){
    miclose_volume(*volume_iterator);
  }

  if(use_mask){
    miclose_volume(mask_handle);
  }
  
  results.erase(results.begin() + added_results, results.end());
  result_inds.erase(result_inds.begin() + added_results, result_inds.end());
  
  return List::create(_["vals"] = results,
                      _["inds"] = result_inds);
}

