AC_INIT([RMINC], 1.3.0.0)
AC_PREREQ([2.68])
AC_CONFIG_SRCDIR([src/slice_loop_functions.c])
AC_PROG_CC
AC_PROG_CXX
m4_include([inst/tools/ax_lib_hdf5.m4])
m4_include([inst/tools/ch_minc_depends.m4])

dnl This adds a --with-build-path option, adding the directories'
dnl include to CPPFLAGS -I, lib to LDFLAGS and the rpath
AC_ARG_WITH([build-path],
	AS_HELP_STRING([--with-build-path=DIR], [build using DIR/include and DIR/lib]),
	[
		#Build-path code taken from Steve Robbins m4 macro
    		for d in `echo $withval | tr : ' '`; do
       		    test -d $d || AC_MSG_ERROR([build path $d not found.])
       		    test -d $d/include && CPPFLAGS="$CPPFLAGS -I$d/include/"
       		    test -d $d/lib && LDFLAGS="$LDFLAGS -L$d/lib/ -Wl,-rpath,$d/lib/"
   		done	
	])


dnl The adds a special variable MINC_TOOLKIT_BUILD_DIR to the configure script
AC_ARG_VAR([MINC_TOOLKIT_BUILD_DIR],[where to build minctoolkit, default /opt/minc-itk4/])

dnl Adds a --enable-build-minc toggle to the configure script
dnl if the switch is not toggled check for minc's dependencies
dnl otherwise check if minc2 is found, if so add it's dependencies 
dnl if not build minc toolkit and try to add the dependencies afterward
AC_ARG_ENABLE(
	[build-minc],
	AS_HELP_STRING([--disable-build-minc], [enable RMINC to build minc-toolkit if it is missing])
)

TOOLKIT_SEARCH

AS_IF([test "x$enable_build_minc" != "xno"], [
	    AC_CHECK_LIB(minc2, mifree_names, 
	         [CHECK_MINC_DEPENDS],
	         [INSTALL_MINC])],
	    [CHECK_MINC_DEPENDS]
) 

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC

dnl Now find the compiler and compiler flags to use
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi

dnl substitute CPPFLAGS and LIBS
AC_SUBST(CPPFLAGS)
AC_SUBST(LIBS)
AC_SUBST(LDFLAGS)
dnl and do subsitution in the src/Makevars.in
AC_OUTPUT(src/Makevars)
