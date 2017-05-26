#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]


// Weighted Quick Union Class, stores a vector of parent indices and subtree sizes
// The weighting decreases tree height by binding smaller trees to larger trees
// Uses path flattening to ensure trees stay flat and root search stays fast
// https://www.cs.princeton.edu/~rs/AlgsDS07/01UnionFind.pdf
class WeightedQuickUnion
{
private:
  std::vector<int> id; 
  std::vector<int> sz;
  std::vector<double> mass;
  int root(int p);
  
public:
  WeightedQuickUnion(int N, std::vector<double> mass);
  bool find(int p, int q);
  void unite(int p, int q);
  int get_root(int p);
  std::vector<int> sizes();
  int sizes(int p);
  double cluster_mass(int p);
};

WeightedQuickUnion::WeightedQuickUnion(int N, std::vector<double> M){
  id = std::vector<int>(N);
  sz = std::vector<int>(N,1);
  mass = M;
  std::iota(id.begin(), id.end(), 0);
}

// Find the root of node i, flattening the path as you go
int WeightedQuickUnion::root(int i){
  while(i != id[i]){
    id[i] = id[id[i]]; // path flattening
    i = id[i];
  }
  return i;
}

// Check if roots are identical, not used for TFCE
bool WeightedQuickUnion::find(int p, int q){
  return root(p) == root(q);
}

// Merge two trees if they don't yet share a root
void WeightedQuickUnion::unite(int p, int q){
  int i = root(p);
  int j = root(q);
  
  if(i != j){
    if(sz[i] < sz[j]){
      id[i] = j; sz[j] += sz[i]; mass[j] += mass[i];
    } else {
      id[j] = i; sz[i] += sz[j]; mass[i] += mass[j];
    }
  }
}

// Find the root of node p without path flattening, non-private
int WeightedQuickUnion::get_root(int p){ 
  while(p != id[p]){
    p = id[p];
  } 
  return p;
}

// Get the size of the subtree rooted at node p
int WeightedQuickUnion::sizes(int p){ return sz[p]; }
double WeightedQuickUnion::cluster_mass(int p){ return mass[p]; }

//Sorter, thanks stack overflow
template <typename T>
std::vector<int> sort_indexes_desc(const std::vector<T> &v) {
  
  // initialize original index locations
  std::vector<int> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);
  
  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end()
         , [&v](int i1, int i2){return v[i1] > v[i2]; });
  
  return idx;
}

// One sided tfce on arbitrary adjacency lists
// [[Rcpp::export]]
std::vector<double> graph_tfce_wqu(std::vector<double> map, std::vector<std::vector<int> > adjacencies
                                  , double E, double H, int nsteps
                                  , std::vector<double> weights){
  //, double nsteps){
  if(map.size() != adjacencies.size()) stop("Mismatch between adjacency list and data");
  
  int nverts = map.size(); 
  std::vector<int> indices = sort_indexes_desc(map);
  //WeightedQuickUnion sets(indices);
  WeightedQuickUnion sets(nverts, weights);
  
  double hmin = 0; //map[indices.back()];
  double hmax = map[indices[0]];
  double dh = (hmax - hmin)/ static_cast<double>(nsteps); //stepsize;
  std::vector<double> tfce = std::vector<double>(nverts, 0);
  if(hmax <= hmin) return tfce;
  
  for(double t = hmax; t >= hmin; t -= dh){
    //Optional threshold tracker
    //Rcout << t << "\n";
    checkUserInterrupt();
    //std::list<int> visited;
    int nvisited = 0;
    for(int sort_ord = 0; sort_ord < nverts && map[indices[sort_ord]] >= t; ++sort_ord){
      int ind = indices[sort_ord];
      // visited.push_back(ind); //.insert(ind);
      ++nvisited;
      for(int neb = 0; neb < adjacencies[ind].size(); ++neb){
        int neighbour = adjacencies[ind][neb]; 
        double neighbour_val = map[neighbour];
        if(neighbour_val >= t && neighbour_val <= map[ind]){//&& visited.count(neighbour) == 0){
          //Find is no longer necessary because unite checks if the roots are identical
          //if(!sets.find(neighbour, ind)){
            sets.unite(neighbour, ind);
          //}
        }
      }
    }
    
    double height_mul = pow(t,H);
    //for(std::list<int>::iterator it = visited.begin(); it != visited.end(); ++it){
    for(int i = 0; i <= nvisited; ++i){
      //int ind = *it;
      int ind = indices[i];
      int root = sets.get_root(ind);
      double tree_size = sets.cluster_mass(root); //mass of root nodes are cluster sizes
      tfce[ind] += (height_mul * pow(tree_size, E) * dh);
    }
    
    
  }
  
  return tfce;
}

// Two sided tfce on arbitrary adjacency graphs
// [[Rcpp::export]]
std::vector<double> graph_tfce(std::vector<double> map, std::vector<std::vector<int> > adjacencies
                                     , double E, double H, int nsteps
                                     , std::vector<double> weights){
 std::vector<double> pos = graph_tfce_wqu(map, adjacencies, E, H, nsteps, weights);
 std::transform(map.begin(), map.end(), map.begin(), std::bind1st(std::multiplies<double>(), -1));
 std::vector<double> neg = graph_tfce_wqu(map, adjacencies, E, H, nsteps, weights);
 for(int i = 0; i < pos.size(); ++i) pos[i] -= neg[i];
 return pos;
}

//////////////////////////////////////////////////////////
// Code for generating adjacency graphs for regular grids


// [[Rcpp::export]]
int coords2ind(int i, int j, int k, int d1, int d2, int d3){
  return k * d2 * d1 + j * d1 + i; 
}

// [[Rcpp::export]]
std::vector<int> ind2coords(int v, int d1, int d2, int d3){
  std::vector<int> res(3);
  res[2] = v / (d1 * d2);
  res[1] = v % (d1 * d2) / d1;
  res[0] = v % d1;
  return res;
} 

// [[Rcpp::export]]
std::vector<std::vector<int> > neighbour_list(double x, double y, double z, int n) {
  int maxdist;
  
  // Set maximum distance a neighbour can be
  switch(n){
    case 6: maxdist = 1; break;
    case 18: maxdist = 2; break;
    case 26: maxdist = 3; break;
    default: stop("Nearest neighbour connectivity must be 6,18, or 26");
  }
  
  int nv = x * y * z; 
  std::vector<std::vector<int> > nlist(nv);
  
  for(int i = 0; i < x; ++i){
    checkUserInterrupt();
    for(int j = 0; j < y; ++j){
      for(int k = 0; k < z; ++k){
        std::vector<int> nebs;
        nebs.reserve(n);
        for(int oi = -1; oi < 2; ++oi){
          if(i + oi >= 0 && i + oi < x){
            for(int oj = -1; oj < 2; ++oj){
              if(j + oj >= 0 && j + oj < y){
                for(int ok = -1; ok < 2; ++ok){
                  if(k + ok >= 0 && k + ok < z){
                    if(abs(oi) + abs(oj) + abs(ok) <= maxdist && abs(oi) + abs(oj) + abs(ok) != 0){
                      nebs.push_back(coords2ind(i + oi, j + oj, k + ok, x,y,z));
                    }
                  }
                }
              }
            }
          }
        }
        nlist[coords2ind(i, j, k, x,y,z)] = nebs;
      }
    }
  }
  
  return nlist;
}

// [[Rcpp::export]]
std::vector<double> mesh_area(std::vector<double> vertices, std::vector<double> triangles){
  std::vector<double> areas(vertices.size()/3, 0);
  
  int linear_ind = 0;
  int v1_ind, v2_ind, v3_ind;
  double d1x, d1y, d1z;
  double d2x, d2y, d2z;
  double area;

  for(int i = 0; i < (triangles.size() / 3); ++i){
    linear_ind = 3 * i;
    v1_ind = 3 * (triangles[linear_ind] - 1);
    v2_ind = 3 * (triangles[linear_ind + 1] - 1);
    v3_ind = 3 * (triangles[linear_ind + 2] - 1);
    
    d1x = vertices[v2_ind] - vertices[v1_ind];
    d1y = vertices[v2_ind + 1] - vertices[v1_ind + 1];
    d1z = vertices[v2_ind + 2] - vertices[v1_ind + 2];
    
    d2x = vertices[v3_ind] - vertices[v1_ind];
    d2y = vertices[v3_ind + 1] - vertices[v1_ind + 1];
    d2z = vertices[v3_ind + 2] - vertices[v1_ind + 2];
    
    area = .5 * sqrt(
      pow(d1y*d2z - d1z*d2x, 2) +
        pow(d1z*d2x - d1x*d2z, 2) +
        pow(d1x*d2y - d1y*d2x, 2));
    
    areas[v1_ind/3] += (1/3.0) * area;
    areas[v2_ind/3] += (1/3.0) * area;
    areas[v3_ind/3] += (1/3.0) * area;
  }
  
  return areas;
}

