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
               AS_IF([test x$HOME != "x"],
                     AS_IF([test -d $HOME/local/minc-itk4/], 
                           [TOOLKIT_PATH="$HOME/local/minc-itk4/"]))
               AS_IF([test x$MINC_TOOLKIT_BUILD_PATH != "x"],
                     AS_IF([test -d $MINC_TOOLKIT_BUILD_PATH],  
                           [TOOLKIT_PATH="$MINC_TOOLKIT_BUILD_PATH"]))
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
        AX_LIB_HDF5
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


     
dnl	AS_IF([test -z $with_build_path && test x$enable_build_minc != "xyes"],
dnl              [
dnl                 TOOLKIT_SEARCH
dnl                 AC_CHECK_HEADER(hdf5.h, , [
dnl                     AC_CHECK_PROG(H5CC, h5cc, "yes", "no")
dnl	             AS_IF(test x$H5CC = "xyes", [
dnl	                AC_REQUIRE([AC_PROG_GREP])
dnl	                H5CC_LDFLAGS=$(h5cc -show | eval "$GREP -o '\-L[[^ ]]*'" | tr "\n" " ") 
dnl	                H5CC_CPPFLAGS=$(h5cc -show | eval "$GREP -o '\-I[[^ ]]*'" | tr "\n" " ")
dnl	                LDFLAGS="$LDFLAGS $H5CC_LDFLAGS"
dnl	                CPPFLAGS="$CPPFLAGS $H5CC_CPPFLAGS"
dnl	              ])
dnl                  ])	  
dnl              ])
dnl               
