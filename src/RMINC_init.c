#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void convert_voxel_to_world(void *, void *, void *, void *, void *);
extern void convert_world_to_voxel(void *, void *, void *, void *, void *);
extern void get_hyperslab(void *, void *, void *, void *);
extern void get_volume_sizes(void *, void *);
extern void get_voxel_from_files(void *, void *, void *, void *, void *, void *);
extern void get_world_voxel_from_files(void *, void *, void *, void *, void *, void *);
extern void qvalue_min(void *, void *, void *, void *);
extern void write_minc2_volume(void *, void *, void *, void *, void *, void *, void *);

/* .Call calls */
extern SEXP get_minc_history(SEXP);
extern SEXP get_minc_separations(SEXP, SEXP);
extern SEXP get_vector_from_files(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP minc2_model(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP minc_overwrite_history(SEXP, SEXP, SEXP);
extern SEXP per_voxel_anova(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _RMINC_anat_summary(SEXP, SEXP, SEXP);
extern SEXP _RMINC_coords2ind(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _RMINC_count_labels(SEXP);
extern SEXP _RMINC_graph_tfce(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _RMINC_graph_tfce_wqu(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _RMINC_ind2coords(SEXP, SEXP, SEXP, SEXP);
extern SEXP _RMINC_mesh_area(SEXP, SEXP);
extern SEXP _RMINC_neighbour_list(SEXP, SEXP, SEXP, SEXP);
extern SEXP _RMINC_rcpp_minc_apply(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _RMINC_replaceValues(SEXP, SEXP, SEXP, SEXP);
extern SEXP vertex_anova_loop(SEXP, SEXP, SEXP);
extern SEXP vertex_lm_loop(SEXP, SEXP, SEXP);
extern SEXP vertex_wlm_loop(SEXP, SEXP, SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
  {"convert_voxel_to_world",     (DL_FUNC) &convert_voxel_to_world,     5},
  {"convert_world_to_voxel",     (DL_FUNC) &convert_world_to_voxel,     5},
  {"get_hyperslab",              (DL_FUNC) &get_hyperslab,              4},
  {"get_volume_sizes",           (DL_FUNC) &get_volume_sizes,           2},
  {"get_voxel_from_files",       (DL_FUNC) &get_voxel_from_files,       6},
  {"get_world_voxel_from_files", (DL_FUNC) &get_world_voxel_from_files, 6},
  {"qvalue_min",                 (DL_FUNC) &qvalue_min,                 4},
  {"write_minc2_volume",         (DL_FUNC) &write_minc2_volume,         7},
  {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
  {"get_minc_history",       (DL_FUNC) &get_minc_history,        1},
  {"get_minc_separations",   (DL_FUNC) &get_minc_separations,    2},
  {"get_vector_from_files",  (DL_FUNC) &get_vector_from_files,   6},
  {"minc2_model",            (DL_FUNC) &minc2_model,            11},
  {"minc_overwrite_history", (DL_FUNC) &minc_overwrite_history,  3},
  {"per_voxel_anova",        (DL_FUNC) &per_voxel_anova,         7},
  {"_RMINC_anat_summary",     (DL_FUNC) &_RMINC_anat_summary,      3},
  {"_RMINC_coords2ind",       (DL_FUNC) &_RMINC_coords2ind,        6},
  {"_RMINC_count_labels",     (DL_FUNC) &_RMINC_count_labels,      1},
  {"_RMINC_graph_tfce",       (DL_FUNC) &_RMINC_graph_tfce,        6},
  {"_RMINC_graph_tfce_wqu",   (DL_FUNC) &_RMINC_graph_tfce_wqu,    6},
  {"_RMINC_ind2coords",       (DL_FUNC) &_RMINC_ind2coords,        4},
  {"_RMINC_mesh_area",        (DL_FUNC) &_RMINC_mesh_area,         2},
  {"_RMINC_neighbour_list",   (DL_FUNC) &_RMINC_neighbour_list,    4},
  {"_RMINC_rcpp_minc_apply",  (DL_FUNC) &_RMINC_rcpp_minc_apply,  10},
  {"_RMINC_replaceValues",    (DL_FUNC) &_RMINC_replaceValues,     4},
  {"vertex_anova_loop",      (DL_FUNC) &vertex_anova_loop,       3},
  {"vertex_lm_loop",         (DL_FUNC) &vertex_lm_loop,          3},
  {"vertex_wlm_loop",        (DL_FUNC) &vertex_wlm_loop,         4},
  {NULL, NULL, 0}
};

void R_init_RMINC(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
