#include <Rcpp.h>
#include "minc2.h"
#include "minc_cpp.h"
#include "minc_reader.h"
#include <sstream>
#include <stdlib.h>
#include <math.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
List anat_summary(CharacterVector filenames,
                  IntegerVector atlas,
                  std::string method,
                  LogicalVector mask,
                  bool use_mask) {

  vector<mihandle_t> volumes = open_minc2_volumes(filenames);

  if(!check_same_dimensions(volumes)){
    for(int i = 0; i < volumes.size(); ++i){
      miclose_volume(volumes[i]);
    }
    stop("At least one file is a different size");
  }

  vector<misize_t> sizes = get_volume_dimensions(volumes[0]);
  int total_voxels = sizes[0] * sizes[1] * sizes[2];

  misize_t vol_sizes[3];
  for(int i = 0; i < 3; ++i){
    vol_sizes[i] = sizes[i];
  }
  misize_t offsets[3] = {0,0,0};

  if(atlas.size() != total_voxels){
    stop("Atlas and files differ in size");
  }
  
  if(use_mask && mask.size() != total_voxels){
    stop("Mask and files differ in size");
  }

  int max_label = 0;
  for(int at_vox = 0; at_vox < total_voxels; ++at_vox){
    if(atlas[at_vox] > max_label) max_label = atlas[at_vox];
  }

  IntegerMatrix label_counts(max_label + 1, filenames.size());
  NumericMatrix label_values(max_label + 1, filenames.size());

  double *vol_buffer =
    (double *) malloc(total_voxels * sizeof(double));
  
  NumericVector separations = as<NumericVector>(get_minc_separations(wrap(filenames[0])));
  double vox_vol = fabs(separations[0]) * fabs(separations[1]) * fabs(separations[2]);

  for(int subject = 0; subject < volumes.size(); ++subject){
    stringstream error_message;
    error_message << "Error Reading Volume " << (subject + 1) << "\n";

    cautious_get_hyperslab(volumes[subject],
                           MI_TYPE_DOUBLE,
                           offsets,
                           vol_sizes,
                           vol_buffer,
                           error_message.str());

    for(int i = 0; i < total_voxels; ++i){
      int current_label = atlas[i];
      double current_value = vol_buffer[i];
      
      if(!use_mask || mask[i]){
        ++label_counts(current_label, subject);
        if(method == "jacobians"){
          label_values(current_label, subject) += exp(current_value) * vox_vol;  
        } else {
          label_values(current_label, subject) += current_value;
        }
      }
    }
  }
  
  free(vol_buffer);
  for(int i = 0; i < volumes.size(); ++i){
    miclose_volume(volumes[i]);
  }
  
  return(List::create(_["values"] = label_values, _["counts"] = label_counts));
}


// [[Rcpp::export]]
NumericMatrix count_labels(CharacterVector filenames,
                           IntegerVector atlas,
                           LogicalVector mask,
                           bool use_mask) {
  
  vector<mihandle_t> volumes = open_minc2_volumes(filenames);
  
  if(!check_same_dimensions(volumes)){
    for(int i = 0; i < volumes.size(); ++i){
      miclose_volume(volumes[i]);
    }
    stop("At least one file is a different size");
  }
  
  vector<misize_t> sizes = get_volume_dimensions(volumes[0]);
  int total_voxels = sizes[0] * sizes[1] * sizes[2];
  
  misize_t vol_sizes[3];
  for(int i = 0; i < 3; ++i){
    vol_sizes[i] = sizes[i];
  }
  misize_t offsets[3] = {0,0,0};
  
  if(atlas.size() != total_voxels){
    stop("Atlas and files differ in size");
  }
  
  if(use_mask && mask.size() != total_voxels){
    stop("Mask and files differ in size");
  }
  
  int max_label = 0;
  for(int at_vox = 0; at_vox < total_voxels; ++at_vox){
    if(atlas[at_vox] > max_label) max_label = atlas[at_vox];
  }
  
  NumericMatrix label_values(max_label + 1, filenames.size());
  
  double *vol_buffer =
    (double *) malloc(total_voxels * sizeof(double));
  
  NumericVector separations = as<NumericVector>(get_minc_separations(wrap(filenames[0])));
  double vox_vol = fabs(separations[0]) * fabs(separations[1]) * fabs(separations[2]);
  
  for(int subject = 0; subject < volumes.size(); ++subject){
    stringstream error_message;
    error_message << "Error Reading Volume " << (subject + 1) << "\n";
    
    cautious_get_hyperslab(volumes[subject],
                           MI_TYPE_DOUBLE,
                           offsets,
                           vol_sizes,
                           vol_buffer,
                           error_message.str());
    
    for(int i = 0; i < total_voxels; ++i){
      int current_label = (int)(vol_buffer[i] + 0.5);
      
      if(current_label > label_values.nrow()){
        stop("A label not present in the atlas was found in a volume");
      }
      
      if(!use_mask || mask[i]){
        label_values(current_label, subject) += vox_vol;  
      } 
    }
  }
  
  free(vol_buffer);
  for(int i = 0; i < volumes.size(); ++i){
    miclose_volume(volumes[i]);
  }
  
  return(label_values);
}  
