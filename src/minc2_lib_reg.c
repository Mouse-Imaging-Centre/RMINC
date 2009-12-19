/*
	Update: Increasing the reference count on the HDF5 library caused
			aborts.  So, I decided to:
			(1) dlopen the minc2 lib, since any this should result in an
				increment of the hdf5 lib (since minc2 links to it)
				
			(2) do *not* decrement on de-registration.  Apparently, we *need*
				the additional reference. Dunno why, but it works.
				
*/


#include "minc2_io.h"
#include <R_ext/Rdynload.h>
#include <dlfcn.h>								// need this for dlopen, dlclose, etc



/*	define a static pointer to the hdf5 reference, since the registration
	routine needs to pass it to the de-registration routine.   */
static void		*hdf_dl_chkptr = 0;

void R_init_RMINC(DllInfo *info) {

//	void		*hdf_dl_ptr;


	/* register routine, and allocate resources */
	if ( R_DEBUG_mincIO ) Rprintf("Registering shared library RMINC in R_init_RMINC\n");


	// load the hdf5 shared library (just to see if we can)
	if ( R_DEBUG_mincIO ) Rprintf("Pointer to library prior to dlopen call: %lu\n", hdf_dl_chkptr);
	if ( R_DEBUG_mincIO ) Rprintf("Attempting to load the minc2 library ... ");
	hdf_dl_chkptr = dlopen("libminc2.dylib", RTLD_LOCAL | RTLD_LAZY);
	// if ( !hdf_dl_chkptr ) {
	// 	Rprintf ("Nope. \nError: %s\n", dlerror());
	// 	exit(1);
	// } else {
	// 	Rprintf("Loaded!\n");
	// }
	if ( R_DEBUG_mincIO ) Rprintf("Pointer to library: %x\n", hdf_dl_chkptr);

}


void R_unload_RMINC(DllInfo *info) {

	int			hdf_dl_status;


	/* release resources */
	if ( R_DEBUG_mincIO ) Rprintf("Unloading shared library RMINC in R_unload_RMINC\n");


	// first, we need to get the pointer to the library
	// ... so check to see if the hdf5 library is loaded, since
	// ... this will return the pointer if it is (else NULL)
	if ( R_DEBUG_mincIO ) Rprintf("Pointer to library: %lu\n", hdf_dl_chkptr);
	if ( R_DEBUG_mincIO ) Rprintf("Is the minc2 library loaded ? ... ");
	

//	hdf_dl_ptr = dlopen("libhdf5.dylib", RTLD_LOCAL | RTLD_LAZY | RTLD_NOLOAD);
if ( R_DEBUG_mincIO ) {
	if ( !hdf_dl_chkptr ) {
		// no, nothing to do
		Rprintf ("No.  Nothing to unload, I suppose.\n");
		//
	} else {
		// yes, let's unload it (or, perhaps, do nothing)
		Rprintf("Yes!\n");
		// Rprintf("UnLoading the hdf library ... ");
		// hdf_dl_status = dlclose(hdf_dl_chkptr);
		// if ( hdf_dl_status ) { 
		// 	Rprintf("\nError unloading hdf5 library. Status %d\n", hdf_dl_status);
		// //		exit(1);
		// } else {
		// printf("Done.\n");
		// };
	}
}


}



