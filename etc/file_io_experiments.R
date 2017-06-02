library(inline)

local_cpp_flags <-
  " -I/include -I/axiom2/projects/software/arch/linux-precise/include"

local_lib_flags <-
  "-L/projects/mice/share/arch/linux-3_2_0-36-generic-x86_64-eglibc-2_15/lib  -L/axiom2/projects/software/arch/linux-precise/lib -Wl,-rpath,/axiom2/projects/software/arch/linux-precise/lib ${LAPACK_LIBS} ${BLAS_LIBS} -lminc2  -lm -lrt -lz  -lhdf5 -lhdf5_hl ${FLIBS}"

hpf_cpp_flags <- 
  "-I/hpf/largeprojects/MICe/tools/minc-toolkit/1.9.11/include"

hpf_lib_flags <-
  "-L/hpf/largeprojects/MICe/tools/minc-toolkit/1.9.11/lib -Wl,-rpath,/hpf/largeprojects/MICe/tools/minc-toolkit/1.9.11/lib ${LAPACK_LIBS} ${BLAS_LIBS} -lminc2   ${FLIBS}"

src<-'
mihandle_t hvol;
misize_t start[3] = {0,0,0};
misize_t count[3] = {9,9,9};

double* buffer = (double *) malloc(sizeof(double) * 1000); 

char* file = "/hpf/largeprojects/MICe/tools/rminctestdata//testsummaryminc10.mnc";
int res = miopen_volume(file, MI2_OPEN_READ, &hvol);

if(res != MI_NOERROR) error("read failed");

res = miget_real_value_hyperslab(hvol, MI_TYPE_DOUBLE, start, count, buffer);

for(int i = 0; i < 1000; ++i)
  x[i] = buffer[i];

free(buffer);
miclose_volume(hvol);
'

big_read <-
  cfunction(sig = c(x = "numeric")
            , body = src
            , cppargs = hpf_cpp_flags
            , libargs = hpf_lib_flags
            , includes = '#include <minc2.h>'
            , convention = ".C"
            , language = "C")

x_test <- numeric(1000)

big_read(x_test)

src2 <-'
mihandle_t hvol;
misize_t start[3];
misize_t count[3] = {0,0,1};

double* buffer = (double *) malloc(sizeof(double) * 1); 

char* file = "/hpf/largeprojects/MICe/tools/rminctestdata//testsummaryminc10.mnc";
int res = miopen_volume(file, MI2_OPEN_READ, &hvol);

if(res != MI_NOERROR) error("read failed");

int i = 0;
for(misize_t x = 0; x < 10; ++x)
  for(misize_t y = 0; y < 10; ++y)
    for(misize_t z = 0; z < 9; ++z){
      start[0] = x; start[1] = y; start[2] = z;
      res = miget_real_value_hyperslab(hvol, MI_TYPE_DOUBLE, start, count, buffer);
      v[i] = buffer[0];
      ++i;
    }

free(buffer);
miclose_volume(hvol);
'

v_test <- numeric(1000)

small_reads <-
  cfunction(sig = c(v = "numeric"),
            , body = src2
            , cppargs = hpf_cpp_flags
            , libargs = hpf_lib_flags
            , includes = '#include "minc2.h"'
            , convention = ".C"
            , language = "C")

small_reads(v_test)

microbenchmark::microbenchmark(
  small_reads(v_test)
  , big_read(x_test)
)
