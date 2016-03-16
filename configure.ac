AC_INIT([RMINC], 1.3.0.0)
AC_PREREQ([2.68])
AC_CONFIG_SRCDIR([src/slice_loop_functions.c])


dnl This adds a --with-build-path option, adding the directories'
dnl include to CPPFLAGS -I, lib to LDFLAGS and the rpath
AC_ARG_WITH([build-path],
	AS_HELP_STRING([--with-build-path=DIR], [build using DIR/include and DIR/lib]),
	[
		#Build-path code taken from Steve Robbins m4 macro
    		for d in `echo $withval | tr : ' '`; do
       		    test -d $d || AC_MSG_ERROR([build path $d not found.])
       		    test -d $d/include && CPPFLAGS="$CPPFLAGS -I$d/include"
       		    test -d $d/lib && LDFLAGS="$LDFLAGS -L$d/lib -Wl,-rpath=$d/lib"
   		done
	])

dnl The adds a special variable MINC_TOOLKIT_BUILD_DIR to the configure scrip
AC_ARG_VAR([MINC_TOOLKIT_BUILD_DIR],[where to build minctoolkit, default /opt/minc-itk4/])

dnl This checks for the minc dependencies, called after minc2 is located or installed
AC_DEFUN([CHECK_MINC_DEPENDS],
[
	AC_CHECK_HEADER(zlib.h, , [AC_MSG_ERROR([zlib not found])])
	AC_CHECK_LIB(hdf5, H5Fopen, , [AC_MSG_ERROR([hdf5 not found])])
	AC_CHECK_LIB(minc2, mifree_name, , [AC_MSG_ERROR([minc2 not found])]) 
])


dnl This macro allows RMINC to install minc_toolkit,
dnl which includes minc2 along with hdf5 and zlib dependencies
AC_DEFUN([INSTALL_MINC],
  [
	#Install minc_toolkit
	
	if [[ -z $MINC_TOOLKIT_BUILD_DIR ]];
	   then MINC_TOOLKIT_BUILD_DIR=$HOME/local/minc-itk4/
	fi
	
	if [[ -z $TMPDIR ]];
	   then TMPDIR=/tmp
	fi

	orig_dir=$(pwd)
	cd $TMPDIR
	
	if [[ ! -d minc-toolkit-v2 ]]; then
	   echo Downloading minc-toolkit-v2
	   git clone --recursive https://github.com/BIC-MNI/minc-toolkit-v2
	fi

	echo Building minc-toolkit-v2 in $MINC_TOOLKIT_BUILD_DIR

	cd minc-toolkit-v2
	
	if [[ ! -d build ]]; then
	   mkdir build
	else
	   rm -r build/*
	fi
	
	cd build
	cmake .. -DCMAKE_INSTALL_PREFIX:PATH=$MINC_TOOLKIT_BUILD_DIR \
                  -DCMAKE_BUILD_TYPE:STRING=Release \
             	  -DMT_BUILD_LITE:BOOL=ON \
                  -DMT_BUILD_SHARED_LIBS:BOOL=ON \
               	  -DMT_USE_OPENMP:BOOL=ON \
      		  -DCPACK_BINARY_DEB:BOOL=ON \
                  -DCPACK_BINARY_NSIS:BOOL=OFF \
               	  -DCPACK_BINARY_RPM:BOOL=OFF \
               	  -DCPACK_BINARY_STGZ:BOOL=OFF \
               	  -DCPACK_BINARY_TBZ2:BOOL=OFF \
               	  -DCPACK_BINARY_TGZ:BOOL=OFF \
               	  -DCPACK_BINARY_TZ:BOOL=OFF \
               	  -DCPACK_SOURCE_TBZ2:BOOL=OFF \
               	  -DCPACK_SOURCE_TGZ:BOOL=OFF \
               	  -DCPACK_SOURCE_TZ:BOOL=OFF \
               	  -DCPACK_SOURCE_ZIP:BOOL=OFF \
               	  -DCPACK_SOURCE_TXZ:BOOL=OFF 
	make all
	make install
	
	cd $orig_dir

        CPPFLAGS="$CPPFLAGS -I$MINC_TOOLKIT_BUILD_DIR/include"
       	LDFLAGS="$LDFLAGS -L$MINC_TOOLKIT_BUILD_DIR/lib -Wl,-rpath=$MINC_TOOLKIT_BUILD_DIR/lib"
		
	CHECK_MINC_DEPENDS	
  ])

AC_PROG_CXX
AC_PROG_CC

dnl Adds a --enable-build-minc toggle to the configure script
dnl if the switch is not toggled check for minc's dependencies
dnl otherwise check if minc2 is found, if so add it's dependencies 
dnl if not build minc toolkit and try to add the dependencies afterward
AC_ARG_ENABLE(
	[build-minc],
	AS_HELP_STRING([--enable-build-minc], [enable RMINC to build minc-toolkit if it is missing]),
	[AC_CHECK_LIB(minc2, mifree_names, 
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
dnl and do subsitution in the src/Makevars.in
AC_OUTPUT(src/Makevars)
