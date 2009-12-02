AC_DEFUN([mni_REQUIRE_MINC],
[
    AC_ARG_WITH([minc2],
    [  --with-minc2            build using minc 2.0 libraries],
    [
    ])
    mni_REQUIRE_LIB(m,[#include <math.h>],[double x = sqrt(3.);])
    mni_REQUIRE_LIB(netcdf,[#include <netcdf.h>],[int i = ncopen("",0);])
    if test "$with_minc2" = "yes"; then
        mni_REQUIRE_LIB(z,[#include <zlib.h>],[int f = compress2;])
        mni_REQUIRE_LIB(hdf5,[#include <hdf5.h>],[int f = H5Fopen("",0,H5P_DEFAULT);])
        mni_REQUIRE_LIB(minc2,[#include <minc.h>],[int i = miicv_create();])
        AC_DEFINE([MINC2], 1, [Define if should build with MINC 2.0])
    else
        mni_REQUIRE_LIB(minc,[#include <minc.h>],[int i = miicv_create();])
    fi
])



AC_DEFUN([mni_REQUIRE_VOLUMEIO],
[
    AC_REQUIRE([mni_REQUIRE_MINC])
    if test "$with_minc2" = "yes"; then
      mni_REQUIRE_LIB(volume_io2,
  	              [#include <volume_io.h>],
                      [Volume vol; 
	 	      Real voxel = 0;
                      Real x = convert_voxel_to_value(vol,voxel);])
    else
      mni_REQUIRE_LIB(volume_io,
  	              [#include <volume_io.h>],
                      [Volume vol; 
	 	      Real voxel = 0;
                      Real x = convert_voxel_to_value(vol,voxel);])
    fi
])

AC_DEFUN([mni_REQUIRE_BICPL],
[
    AC_REQUIRE([mni_REQUIRE_VOLUMEIO])
    mni_REQUIRE_LIB(bicpl,
		    [#include <bicpl.h>],
		    [File_formats format;
                     int n_obj;
                     object_struct** obj_list;
                     Status s = input_graphics_file("",&format,&n_obj,&obj_list)])
])

AC_DEFUN([mni_REQUIRE_EBTKS],
[
    AC_LANG_PUSH(C++)
    mni_REQUIRE_LIB(EBTKS, [#include <EBTKS/Path.h>],[Path path;])
    AC_LANG_POP
])

AC_DEFUN([mni_REQUIRE_OOBICPL],
[
    AC_REQUIRE([mni_REQUIRE_BICPL])

    # the regular expression C library
    smr_REQUIRED_LIB(pcre, pcre_compile, pcre.h)
    
    AC_LANG_PUSH(C++)


    # the C++ wrapper for the pcre library
    AC_MSG_CHECKING([for -lpcre++])
    LIBS="-lpcre++ $LIBS"
    AC_TRY_LINK([#include <pcre++.h>],
                [using namespace pcrepp; Pcre EE("[a-z]", "i");],,
                [AC_MSG_ERROR(cannot find pcre++ library)])
    AC_MSG_RESULT([yes])

    mni_REQUIRE_LIB(oobicpl,[#include <mniVolume.h>],[mniVolume vol;])
    AC_LANG_POP
])


AC_DEFUN([mni_REQUIRE_BICINVENTOR],
[
    AC_REQUIRE([mni_REQUIRE_BICPL])
    AC_REQUIRE([mni_REQUIRE_OOBICPL])
    AC_REQUIRE([mni_REQUIRE_OPENINVENTOR])

    AC_LANG_PUSH(C++)
    mni_REQUIRE_LIB(bicInventor,
                    [#include <bicInventor.h>],
                    [SoSeparator* root = bic_graphics_file_to_iv("foo.iv");])
    AC_LANG_POP
])


