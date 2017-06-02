Rcpp::sourceCpp("vertex_tfce.cpp")
library(RMINC)

#benchmark_file <- "/hpf/largeprojects/MICe/chammill/POND/reconstructions/minc/all_t1/MR160-088-0019-01_002_T1.mnc"
benchmark_file <- "downsamp.mnc"
#adj <- neighbour_list(192,240,256,6)
adj <- neighbour_list(50,50,50, 6)

# Export statements in c++ code produce wrappers so they are callable from within R
system.time(ftfce <- fast_tfce(mincGetVolume(benchmark_file), adj, E = .5, H = 2, nsteps = 100))
system.time(tfce <- graph_tfce_wqu(mincGetVolume(benchmark_file), adj, E = .5, H = 2, nsteps = 100))
#mtfce <- mincTFCE(mincGetVolume(benchmark_file), side = "positive", d = "10")