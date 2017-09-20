#include <Rcpp.h>
using namespace Rcpp;

// a fast way to replace values based on their labels
// labelVol : the volume containing the labels
// out : the volume whose values will be filled in
// labels : label values to be replaced
// values: replacement values. labels and values must be the same length

// [[Rcpp::export]]
void replaceValues(NumericVector labelVol, 
                   NumericVector out,
                   NumericVector labels,
                   NumericVector values) {
  // build a map
  std::map<int,double> lmap;
  for (int i=0; i<labels.length(); i++) {
    lmap[labels[i]] = values[i];
  }
  
  // substitute at every voxel
  for (int i=0; i<labelVol.length(); i++){
    out[i] = lmap[labelVol[i]];
  }
}


