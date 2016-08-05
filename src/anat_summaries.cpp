// #include <Rcpp.h>
// #include "minc2.h"
// #include "minc_cpp.h"
// #include <sstream>
// #include <stdlib.h>
// #include <math.h>
// using namespace Rcpp;
// using namespace std;
// 
// List anat_summary(CharacterVector filenames,
//                   IntegerVector atlas,
//                   String method) {
//   
//   vector<mihandle_t> volumes = open_minc2_volumes(filenames);
//   
//   if(!check_same_dimensions(volumes)){
//     for(int i = 0; i < volumes.size(); ++i){
//       miclose_volume(volumes[i]);
//     }
//     stop("At least one file is a different size");
//   }
//   
//   vector<misize_t> sizes = get_volume_dimensions(volumes[0]);
//   int total_voxels = sizes[0] * sizes[1] * sizes[2];
//   
//   misize_t vol_sizes[3];
//   for(int i = 0; i < 3; ++i){
//     vol_sizes[i] = sizes[i];
//   }
//   misize_t offsets[3] = {0,0,0};
//   
//   if(atlas.size() != total_voxels){
//     stop("Atlas and files differ in size");
//   }
//   
//   set<int> unique_labels;
//   set<int>::iterator label_it;
//   for(int i = 0; i < atlas.size(); ++i){
//     unique_labels.insert(atlas[i]);
//   }
//   int max_label;
//   for(label_it = unique_labels.begin(); label_it != unique_labels.end(); ++label_it){
//     if(*label_it > max_label) max_label = *label_it;
//   }
//   
//   IntegerVector label_counts(max_label, 0);
//   NumericVector label_values(max_label, 0);
//   NumericMatrix output_label_values(max_label, filenames.size());
// 
//   double *vol_buffer =
//     (double *) malloc(total_voxels * sizeof(double));
//   
//   for(int subject; subject < volumes.size(); ++subject){
//     map<double, double> subject_res;
//     
//     stringstream error_message;
//     error_message << "Error Reading Volume " << (subject + 1) << "\n";
//     
//     cautious_get_hyperslab(volumes[subject],
//                            MI_TYPE_DOUBLE,
//                            offsets,
//                            vol_sizes,
//                            vol_buffer,
//                            error_message.str());
//     
//     for(int i = 0; i < total_voxels; ++i){
//       
//     } 
//   }
//     
// 
//   
//   