dnl @synopsis smr_OPTIONAL_LIB(LIBRARY,FUNCTION,HEADER-FILE,[OTHER-LIBRARIES],[INCLUDES])
dnl Check for presence of library and header file.
dnl If found, the AC_SUBST(LIBRARY,HAVE_[LIBRARY]) is called.
dnl
dnl @version $Id: smr_OPTIONAL_LIB.m4,v 1.2 2001/07/25 04:22:49 stever Exp $
dnl @author Steve M. Robbins <smr@debian.org>


AC_DEFUN([smr_OPTIONAL_LIB],
[
    AC_CHECK_HEADERS([$3],
                     HAVE_lib[$1]=yes,
		     [AC_MSG_WARN(cannot find header; $1 not available)
                      HAVE_lib[$1]=no],
		     [$5])

    if test $HAVE_lib[$1] = yes; then
        AC_CHECK_LIB([$1],[$2],,
		     [AC_MSG_WARN(cannot find library; $1 not available)
                      HAVE_lib[$1]=no],
		     [$4])
    fi

    AC_SUBST(HAVE_lib[$1])
])

