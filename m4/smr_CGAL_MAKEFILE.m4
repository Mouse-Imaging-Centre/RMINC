dnl @synopsis smr_CGAL_MAKEFILE
dnl
dnl This macro adds a "--with-cgal-makefile" option to the configure script,
dnl in order to specify where the CGAL skeletal makefile is located.
dnl
dnl @version $Id: smr_CGAL_MAKEFILE.m4,v 1.2 2002/05/07 01:13:25 stever Exp $
dnl @author Steve M. Robbins <smr@debian.org>

AC_DEFUN([smr_CGAL_MAKEFILE],
[
    AC_ARG_WITH(cgal-makefile,
[  --with-cgal-makefile=PATHNAME 
                          specify location of CGAL makefile],
[
    if test ! -r "$withval"; then
	AC_MSG_ERROR([Cannot read CGAL makefile "$withval"])
    else
	CGAL_MAKEFILE=$withval
    fi
])

if test -z "$CGAL_MAKEFILE"; then
    for cm in cgal_makefile ../cgal_makefile ../../cgal_makefile none; do
	test -r $cm && break
    done
    if test $cm = none; then
	AC_MSG_ERROR(
[Cannot find a CGAL makefile fragment to include.
Use --with-cgal-makefile or create a symbolic link named \"cgal_makefile\"
in the top of the source tree.])
    fi
    CGAL_MAKEFILE='$(top_builddir)/'$cm
fi

AC_SUBST(CGAL_MAKEFILE)
])

