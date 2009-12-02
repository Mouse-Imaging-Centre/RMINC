dnl @synopsis smr_REQUIRED_LIB(LIBRARY,FUNCTION,HEADER-FILE,[OTHER-LIBRARIES],[INCLUDES])
dnl
dnl @version $Id: smr_REQUIRED_LIB.m4,v 1.2 2001/07/25 04:22:49 stever Exp $
dnl @author Steve M. Robbins <smr@debian.org>


AC_DEFUN([smr_REQUIRED_LIB],
[

  dnl Define convenient aliases for the arguments since there are so
  dnl many of them and I keep confusing myself whenever I have to edit
  dnl this macro.
  pushdef([smr_library],     $1)
  pushdef([smr_function],    $2)
  pushdef([smr_header],      $3)
  pushdef([smr_other_libs],  $4)
  pushdef([smr_includes],    $5)


  # Raise error if [smr_header] cannot be found.
  AC_CHECK_HEADER(smr_header,,
	          AC_MSG_ERROR(cannot find required header),
		  smr_includes)

  # Raise error if [smr_library] cannot be found.
  dnl 2001-07-24 - smr
  dnl If I quote the first argument ([smr_library]), it is passed verbatim
  dnl into the configure script, unlike all the other quoted smr_* names.
  dnl Where is the bug?
  AC_CHECK_LIB(smr_library,smr_function,,
	       AC_MSG_ERROR(cannot find required library),
	       smr_other_libs)

  popdef([smr_library])
  popdef([smr_function])
  popdef([smr_header])
  popdef([smr_other_libs])
  popdef([smr_includes])
])


