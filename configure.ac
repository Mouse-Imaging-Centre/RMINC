AC_INIT([RMINC], 1.3.0.0)
AC_PREREQ([2.68])
AC_CONFIG_SRCDIR([src/slice_loop_functions.c])
AC_PROG_CC
AC_PROG_CXX

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

AC_DEFUN([TOOLKIT_SEARCH], [
        AC_MSG_NOTICE([Searching for minc toolkit])
	AC_MSG_CHECKING([Can we find a minc tool (mincinfo) on the search path]) 
        MINC_INFO_PATH=$(which mincinfo)
	AS_IF([test x$MINC_INFO_PATH == "x"], [
               AC_MSG_RESULT(no)], [
               TOOLKIT_FOUND="yes"
               TOOLKIT_PATH=${MINC_INFO_PATH%/bin*}
               AC_MSG_RESULT(yes)
             ])

        AS_IF([test x$TOOLKIT_FOUND != "xyes"], [
               AC_MSG_CHECKING([Can we find minc toolkit in common locations])
               AS_IF([test -d /opt/minc/], [TOOLKIT_PATH="/opt/minc"])
               AS_IF([test -d /opt/minc-itk4/], [TOOLKIT_PATH="/opt/minc-itk4/"])
               AS_IF([test x$HOME != "x" -a -d $HOME/local/minc-itk4/], 
                     [TOOLKIT_PATH="$HOME/local/minc-itk4/"])
               AS_IF([test x$MINC_TOOLKIT_BUILD_PATH != "x" -a -d $MINC_TOOLKIT_BUILD_PATH], 
                     [TOOLKIT_PATH="$MINC_TOOLKIT_BUILD_PATH"])
               AS_IF([test x$TOOLKIT_PATH != "x"], [
                      TOOLKIT_FOUND="yes"
                      AC_MSG_RESULT(yes)
                      AC_MSG_NOTICE([found toolkit in $TOOLKIT_PATH])
                      AC_MSG_NOTICE([override by setting the environment variable MINC_TOOLKIT_BUILD_PATH])
                      AC_MSG_NOTICE([or by setting --with-build-path as a configure argument])
                      ], [ 
                      AC_MSG_RESULT(no) 
                     ])
              ])

	  AS_IF([test x$TOOLKIT_FOUND == "xyes"], [
                  LDFLAGS="$LDFLAGS -L${TOOLKIT_PATH}/lib -Wl,-rpath,${TOOLKIT_PATH}/lib"
                  CPPFLAGS="$CPPFLAGS -I${TOOLKIT_PATH}/include"
                ])
         ])

dnl This checks for the minc dependencies, called after minc2 is located or installed
AC_DEFUN([CHECK_MINC_DEPENDS],
[
	AC_CHECK_HEADER(zlib.h, , [AC_MSG_ERROR([zlib not found])])	
     
	AS_IF([test -z $with_build_path && test x$enable_build_minc != "xyes"],
              [
                 TOOLKIT_SEARCH
                 AC_CHECK_HEADER(hdf5.h, , [
                     AC_CHECK_PROG(H5CC, h5cc, "yes", "no")
	             AS_IF(test x$H5CC = "xyes", [
	                AC_REQUIRE([AC_PROG_GREP])
	                H5CC_LDFLAGS=$(h5cc -show | eval "$GREP -o '\-L[[^ ]]*'" | tr "\n" " ") 
	                H5CC_CPPFLAGS=$(h5cc -show | eval "$GREP -o '\-I[[^ ]]*'" | tr "\n" " ")
	                LDFLAGS="$LDFLAGS $H5CC_LDFLAGS"
	                CPPFLAGS="$CPPFLAGS $H5CC_CPPFLAGS"
	              ])
                  ])	  
              ])
               
	AC_CHECK_LIB(hdf5, H5open, , [AC_MSG_ERROR([HDF5 not found])])
        AC_CHECK_LIB(minc2, mifree_name, , [AC_MSG_ERROR([minc2 not found])])
])


dnl This macro allows RMINC to install minc_toolkit,
dnl which includes minc2 along with hdf5 and zlib dependencies
AC_DEFUN([INSTALL_MINC],
  [
	#Install minc_toolkit
	
	AS_IF([test x$MINC_TOOLKIT_BUILD_DIR = "x"], [ 
	   MINC_TOOLKIT_BUILD_DIR=$HOME/local/minc-itk4/
	])
	
	AS_IF([test x$TMPDIR = "x"], [ 
	   TMPDIR=/tmp
	])

	orig_dir=$(pwd)
	cd $TMPDIR
	
	AS_IF([test ! -d minc-toolkit-v2], [
	   echo Downloading minc-toolkit-v2
	   git clone --recursive https://github.com/BIC-MNI/minc-toolkit-v2
	])

	AS_ECHO([Building minc-toolkit-v2 in $MINC_TOOLKIT_BUILD_DIR])

	cd minc-toolkit-v2
	
	AS_IF([test ! -d build], [
	   mkdir build
	], [
	   rm -r build/*
	])
	
	cd build
	cmake ..  -DCMAKE_INSTALL_PREFIX:PATH=$MINC_TOOLKIT_BUILD_DIR \
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
       	LDFLAGS="$LDFLAGS -L$MINC_TOOLKIT_BUILD_DIR/lib -Wl,-rpath,$MINC_TOOLKIT_BUILD_DIR/lib"
		
	CHECK_MINC_DEPENDS	
  ])


dnl Adds a --enable-build-minc toggle to the configure script
dnl if the switch is not toggled check for minc's dependencies
dnl otherwise check if minc2 is found, if so add it's dependencies 
dnl if not build minc toolkit and try to add the dependencies afterward
AC_ARG_ENABLE(
	[build-minc],
	AS_HELP_STRING([--enable-build-minc], [enable RMINC to build minc-toolkit if it is missing])
)

AS_IF([test "x$enable_build_minc" = "xyes"], [
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
