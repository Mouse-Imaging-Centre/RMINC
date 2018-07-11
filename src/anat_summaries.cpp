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
                  std::string method) {

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
  
  int max_label = 0;
  for(int at_vox = 0; at_vox < total_voxels; ++at_vox){
    if(atlas[at_vox] > max_label) max_label = atlas[at_vox];
  }

  IntegerMatrix label_counts(max_label + 1, filenames.size());
  NumericMatrix label_values(max_label + 1, filenames.size());

  double *vol_buffer =
    (double *) malloc(total_voxels * sizeof(double));

  NumericVector vox_vol(filenames.size());

  for(int i = 0; i < vox_vol.size(); ++i){
    NumericVector separations = as<NumericVector>(get_minc_separations(wrap(filenames[i])));
    double sep  = fabs(separations[0]) * fabs(separations[1]) * fabs(separations[2]);
    vox_vol[i] = sep;
  }

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
      ++label_counts(current_label, subject);
      
      if(method == "jacobians"){
        label_values(current_label, subject) += exp(current_value) * vox_vol[subject];  
      } else {
        label_values(current_label, subject) += current_value;
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
NumericMatrix count_labels(CharacterVector filenames) {
  
  vector<mihandle_t> volumes = open_minc2_volumes(filenames);
  
  // if(!check_same_dimensions(volumes)){
  //   for(int i = 0; i < volumes.size(); ++i){
  //     miclose_volume(volumes[i]);
  //   }
  //   stop("At least one file is a different size");
  // }
  
  vector<misize_t> sizes = get_volume_dimensions(volumes[0]);
  int total_voxels = sizes[0] * sizes[1] * sizes[2];
  
  misize_t vol_sizes[3];
  for(int i = 0; i < 3; ++i){
    vol_sizes[i] = sizes[i];
  }
  misize_t offsets[3] = {0,0,0};
  
  double *vol_buffer =
    (double *) malloc(total_voxels * sizeof(double));
  
  NumericVector separations = as<NumericVector>(get_minc_separations(wrap(filenames[0])));
  double vox_vol = fabs(separations[0]) * fabs(separations[1]) * fabs(separations[2]);
  
  // Get first subject
  cautious_get_hyperslab(volumes[0],
                         MI_TYPE_DOUBLE,
                         offsets,
                         vol_sizes,
                         vol_buffer,
                         "trouble reading first volume");
  int max_label = 0;
  for(int at_vox = 0; at_vox < total_voxels; ++at_vox){
    if(vol_buffer[at_vox] > max_label) max_label = (int)(vol_buffer[at_vox] + .5);
  }
  NumericMatrix label_values(max_label + 1, filenames.size());
  
  for(int subject = 0; subject < volumes.size(); ++subject){
    stringstream error_message;
    error_message << "Error Reading Volume " << (subject + 1) << "\n";
    
    sizes = get_volume_dimensions(volumes[subject]);
    total_voxels = sizes[0] * sizes[1] * sizes[2];
    for(int i = 0; i < 3; ++i){
      vol_sizes[i] = sizes[i];
    }
    
    free(vol_buffer);
    vol_buffer = (double *) malloc(total_voxels * sizeof(double));
    
    cautious_get_hyperslab(volumes[subject],
                           MI_TYPE_DOUBLE,
                           offsets,
                           vol_sizes,
                           vol_buffer,
                           error_message.str());
    
    for(int i = 0; i < total_voxels; ++i){
      int current_label = (int)(vol_buffer[i] + 0.5);
      
      if(current_label > label_values.nrow()){
        error_message.str("");
        error_message << "A label (" << current_label << ") was found in a volume " << 
            (subject + 1) << " but not in the first volume\n";
        
        stop(error_message.str());
      }
      
      label_values(current_label, subject) += vox_vol;  
    }
  }
  
  free(vol_buffer);
  for(int i = 0; i < volumes.size(); ++i){
    miclose_volume(volumes[i]);
  }
  
  return(label_values);
}  

// [[Rcpp::export]]
List atlas_get_all(CharacterVector filenames,
                   CharacterVector atlas,
                   std::string method){

  if(method != "labels" && atlas.size() != filenames.size() && atlas.size() > 1){
    stop(string("Number of atlases does not match the number of files. ") +
         "anatGetAll requires either a single atlas or one per subject");
  }

  vector<map<int,double> > all_results;
  MincVolume cur_atlas; 
  shared_ptr<double> at_data;
  shared_ptr<int> int_atlas;
  
  if(atlas.size() != 0){
    cur_atlas = MincVolume(as<string>(atlas[0]));
    at_data = cur_atlas.read_volume(MI_TYPE_DOUBLE);
    int_atlas = shared_ptr<int>(new int[cur_atlas.size()]);
    for(int i = 0; i < cur_atlas.size(); ++i)
      int_atlas.get()[i] = round(at_data.get()[i]);
  }

  for(int i = 0; i < filenames.size(); ++i){
    // If the atlas has changed, update all the necessary stuff
    if(i != 0 && atlas.size() > 1){
      cur_atlas = MincVolume(as<string>(atlas[i]));
      at_data = cur_atlas.read_volume(MI_TYPE_DOUBLE);
      int_atlas = shared_ptr<int>(new int[cur_atlas.size()]);
      for(int i = 0; i < cur_atlas.size(); ++i)
        int_atlas.get()[i] = round(at_data.get()[i]);
    } 

    MincVolume cur_subj(as<string>(filenames[i]));
    double vox_vol = cur_subj.voxel_volume();

    if(atlas.size() != 0 && cur_subj.size() != cur_atlas.size())
      stop("File: " + cur_subj.get_filename() +
           " did not match the size of atlas: " +
           cur_atlas.get_filename());

    if(method != "labels")
      shared_ptr<double> subj_data = cur_subj.read_volume(MI_TYPE_DOUBLE);

    map<int,double> res;
    if(method == "means"){
      res = atlasMeans(cur_subj, int_atlas);
    } else if(method == "sums"){
      res = atlasReduce(cur_subj, int_atlas
                        , [](double acc, double x){ return(acc + x); }
                        );
    } else if(method == "jacobians"){
      res = atlasReduce(cur_subj, int_atlas
                        , [=](double acc, double x){ return(acc + vox_vol * exp(x));  }
                        );
    } else if(method == "labels"){
      res = atlasVol(cur_subj, vox_vol);
    } else {
      stop(string("Unrecognized method, this should have been caught in the ") +
           string("parent R code, please file a bug report."));
    }

    all_results.push_back(res);
  }

  return(mergeCountMaps(all_results));
}

map<int,double> atlasMeans(MincVolume subject, shared_ptr<int> atlas_buf){
  shared_ptr<double> subj_buf = subject.read_volume(MI_TYPE_DOUBLE);
  map <int,pair<double,int>> res;
  
  for(int i = 0; i < subject.size(); ++i){
    int key = atlas_buf.get()[i];
    double val = subj_buf.get()[i];

    pair<double,int>& cur_res = res[key];
    cur_res.first += val;
    cur_res.second += 1;
  }

  map<int,double> full_res;
  for(auto it = res.begin(); it != res.end(); ++it){
    full_res[(*it).first] = (*it).second.first / (*it).second.second;
  }
    
  return(full_res);  
}

map<int,double> atlasReduce(MincVolume subject, shared_ptr<int> atlas_buf
                            , std::function<double(double,double)> fun){
  shared_ptr<double> subj_buf = subject.read_volume(MI_TYPE_DOUBLE);
  map<int,double> res;
  
  for(int i = 0; i < subject.size(); ++i){
    int key = atlas_buf.get()[i];
    double val = subj_buf.get()[i];

    res[key] = fun(res[key], val);
  }
    
  return(res);  
}

map<int,double> atlasVol(MincVolume subject, double vox_vol){
  shared_ptr<double> subj_buf = subject.read_volume(MI_TYPE_DOUBLE);
  map<int,double> res;

  for(int i = 0; i < subject.size(); ++i){
    int key = round(subj_buf.get()[i]);
    res[key] += vox_vol;
  }

  return(res);
}

List mergeCountMaps(vector<map<int,double> > maps){
  set<int> keys;
  map<string,NumericVector> merged;
  
  for(auto vit = maps.begin(); vit != maps.end(); ++vit){
    for(auto mit = (*vit).begin(); mit != (*vit).end(); ++mit){
      keys.insert((*mit).first);
    }
  }

  for(auto vit = maps.begin(); vit != maps.end(); ++vit){
    for(auto kit = keys.begin(); kit != keys.end(); ++kit){
      auto val = (*vit).find(*kit);
      if(val == (*vit).end()){
        merged[to_string(*kit)].push_back(NA_REAL);
      } else {
        merged[to_string(*kit)].push_back((*val).second);
      }        
    }
  }

  return(wrap(merged));
}
