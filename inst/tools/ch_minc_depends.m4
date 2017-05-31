dnl Search for the minc-toolkit
AC_DEFUN([MINC_SEARCH], [
  AC_MSG_NOTICE([Searching for libminc])
  
  AS_IF([test "x$MINC_FOUND" != "xyes"], [
    AS_IF([test x$MINC_PATH != "x"], [
      AC_MSG_NOTICE(MINC_PATH was specified)
      MINC_FOUND="yes"
    ])
  ])

  AS_IF([test "x$MINC_FOUND" != "xyes"], [
    AC_MSG_CHECKING("Checking for $MINC_BUILD_PATH (MINC_BUILD_PATH)")
    AS_IF([test "x$MINC_BUILD_PATH" != "x" && test -d "$MINC_BUILD_PATH"], [
      AC_MSG_RESULT(yes)
      MINC_PATH="$MINC_BUILD_PATH"
      MINC_FOUND="yes"
      ], [ AC_MSG_RESULT(no) ])
  ])

  AS_IF([test "x$MINC_FOUND" != "xyes"], [
    AC_MSG_CHECKING(Checking for /opt/minc-itk4/)
    AS_IF([test -d "/opt/minc-itk4"], [
      AC_MSG_RESULT(yes)
      MINC_FOUND="yes"
      MINC_PATH="/opt/minc-itk4"
      ], [ AC_MSG_RESULT(no) ])
  ])

  AS_IF([test "x$MINC_FOUND" != "xyes"], [
    AC_MSG_CHECKING(Checking for /opt/minc/)
    AS_IF([test -d "/opt/minc/"], [
      AC_MSG_RESULT(yes)  
      MINC_FOUND="yes"
      MINC_PATH="/opt/minc"
      ], [ AC_MSG_RESULT(no) ])
  ])

  AS_IF([test x$MINC_FOUND != "xyes"], [
    AC_MSG_CHECKING([Can we find a minc tool (mincinfo) on the search path]) 
    MINC_INFO_PATH=$(which mincinfo)
    AS_IF([test x$MINC_INFO_PATH == "x"], [
      AC_MSG_RESULT(no)], [
      AC_MSG_RESULT(yes)
      MINC_FOUND="yes"
      MINC_PATH=${MINC_INFO_PATH%/bin*}
    ])
  ])

  AS_IF([test "x$MINC_FOUND" == "xyes"], [
    AC_MSG_NOTICE([MINC_PATH set to $MINC_PATH])
    AC_MSG_NOTICE([override by setting the environment variable MINC_PATH])
    AC_MSG_NOTICE([or by setting --with-build-path as a configure argument])

    LDFLAGS="$LDFLAGS -L${MINC_PATH}/lib -Wl,-rpath,${MINC_PATH}/lib"
    CPPFLAGS="$CPPFLAGS -I${MINC_PATH}/include"
  ])
])


dnl This macro builds libminc if it's not found
AC_DEFUN([INSTALL_LIBMINC], [
    AS_IF([test x$MINC_BUILD_PATH = "x"], [ 
	   MINC_BUILD_PATH=$HOME/local/minc-itk4/
    ])
	
    AS_IF([test x$TMPDIR = "x"], [ 
	   TMPDIR=/tmp
    ])

	orig_dir=$(pwd)
	cd $TMPDIR
	
	AS_IF([test ! -d libminc], [
	   echo Downloading libminc
	   git clone --recursive https://github.com/BIC-MNI/libminc
	])

	AS_ECHO([Building libminc in $MINC_BUILD_PATH])

	cd libminc
	
	
	AS_IF([test ! -d build], [
	   mkdir build
	], [
	   rm -r build/*
	])
	
	cd build

	AC_PROG_SED

	AS_IF([test $(uname) = "Darwin"], SED_INPLACE="$SED -i '' -e", SED_INPLACE="$SED -i")

	$SED_INPLACE  "s|\(\$ENV{HDF5_HOME}/lib\)|\1\n/usr/lib/x86_64-linux-g[]nu/hdf5/serial\n/usr/lib/x86_64-linux-g[]nu/hdf5/parallel/|" \
	              ../cmake-modules/FindHDF5.cmake
	
        $SED_INPLACE  "s|\(\$ENV{HDF5_HOME}/i[]nclude\)|\1\n/usr/i[]nclude/hdf5/serial\n/usr/i[]nclude/parallel/|" \
	              ../cmake-modules/FindHDF5.cmake

	cmake .. -DCMAKE_INSTALL_PREFIX:PATH=$MINC_BUILD_PATH

        make all
	make install
	
	cd $orig_dir

        CPPFLAGS="$CPPFLAGS -I$MINC_BUILD_PATH/include"
       	LDFLAGS="$LDFLAGS -L$MINC_BUILD_PATH/lib -Wl,-rpath,$MINC_BUILD_PATH/lib"
		
])    


dnl This checks for the minc dependencies, called after minc2 is located or installed
AC_DEFUN([CHECK_MINC_DEPENDS],
[
	AC_CHECK_HEADER(zlib.h, , [AC_MSG_ERROR([zlib not found])])
        AX_LIB_HDF5

	CFLAGS="$CFLAGS $HDF5_CFLAGS"
        LDFLAGS="$LDFLAGS $HDF5_LDFLAGS"
        CPPFLAGS="$CPPFLAGS $HDF5_CPPFLAGS"
        LIBS="$LIBS $HDF5_LIBS"

	MINC_SEARCH
	
        AS_IF([test x$MINC_FOUND != "xyes"], [
	   INSTALL_LIBMINC
	])

        AC_CHECK_LIB(minc2, mifree_name, , [AC_MSG_ERROR([minc2 not found])])
])



dnl This macro allows RMINC to install minc_toolkit,
dnl which includes minc2 along with hdf5 and zlib dependencies
AC_DEFUN([INSTALL_MINC_TOOLKIT],
  [
	#Install minc_toolkit
	
        #Check for openmp, set USE_OMP if it is
	AC_OPENMP
        AS_IF([test x$OPENMP_CFLAGS == "x"], USE_OMP="off", USE_OMP="on")
          

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
               	  -DMT_USE_OPENMP:BOOL=$USE_OMP \
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
