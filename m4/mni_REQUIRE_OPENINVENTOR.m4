# Macros for finding OpenInventor and related libraries.

# Qt and SoQt sometimes need functions from X11 libraries.
#
AC_DEFUN([mni_REQUIRE_X11],
[
    AC_REQUIRE([AC_PATH_XTRA])

# Qt requires Xext and X11 (and libm)
#
# SoQt needs the above plus SM, ICE, Xi and Xmu
    if test "$no_x" != yes; then
	CPPFLAGS="$CPPFLAGS $X_CFLAGS"
	LDFLAGS="$LDFLAGS $X_LIBS"
	LIBS="$X_PRE_LIBS -lXmu -lXext -lSM -lICE -lXi -lX11 $X_EXTRA_LIBS $LIBS"
    fi
])    


AC_DEFUN([mni_REQUIRE_GL],
[
    AC_REQUIRE([mni_REQUIRE_X11])
    mni_REQUIRE_LIB(GL,[#include <GL/gl.h>],[glEnd();])
])


AC_DEFUN([mni_REQUIRE_GLU],
[
    AC_REQUIRE([mni_REQUIRE_GL])
    mni_REQUIRE_LIB(GLU,[#include <GL/glu.h>],[gluBeginCurve(0);])
])


AC_DEFUN([mni_REQUIRE_OPENINVENTOR],
[
    AC_REQUIRE([mni_REQUIRE_X11])
    AC_REQUIRE([mni_REQUIRE_GL])

    AC_LANG_PUSH(C++)
    mni_LIBS_save=$LIBS
    LIBS="-lCoin $LIBS"

    AC_MSG_CHECKING([for Open Inventor (Coin)])
    AC_TRY_LINK([#include <Inventor/SoDB.h>],
                [SoDB::init();],
	    	[result=yes],
	    	[result=no])
    AC_MSG_RESULT($result)

    if test "$result" = no; then
    	LIBS="-lInventor $mni_LIBS_save"
    	AC_MSG_CHECKING([for Open Inventor (SGI)])
    	AC_TRY_LINK([#include <Inventor/SoDB.h>],
        	    [SoDB::init();],
	   	    [result=yes],
		    [result=no])
    	AC_MSG_RESULT($result)
    fi
    AC_LANG_POP

    if test "$result" = no; then
    	AC_MSG_ERROR([No Open Inventor library detected.])
    fi
])


AC_DEFUN([mni_REQUIRE_QTGL],
[
    AC_REQUIRE([mni_REQUIRE_X11])
    AC_REQUIRE([mni_REQUIRE_GLU])

    AC_MSG_CHECKING([value of the QTDIR environment variable])
    if test x"$QTDIR" = x""; then
	AC_MSG_RESULT([empty])
        AC_MSG_WARN([QTDIR environment variable not set -- this might be an indication of a problem])
    else
 	AC_MSG_RESULT([$QTDIR])
	LDFLAGS="-L$QTDIR/lib $LDFLAGS"
	CPPFLAGS="-I$QTDIR/include $CPPFLAGS"
    fi
    AC_SUBST([QTDIR], [$QTDIR])
    
    AC_LANG_PUSH(C++)
    mni_REQUIRE_LIB(qt,[#include <qapplication.h>],[QString str;])
    AC_LANG_POP

    AC_PATH_PROGS(MOC, moc, echo, $PATH:$QTDIR/bin)
    AC_PATH_PROGS(UIC, uic, echo, $PATH:$QTDIR/bin)
])
    


dnl This used to have -lXi in it.  Perhaps should AC_REQUIRE the AC_PATH_XTRA
dnl macro and do something intelligent with it.
AC_DEFUN([mni_REQUIRE_SOQT],
[
    AC_REQUIRE([mni_REQUIRE_X11])
    AC_REQUIRE([mni_REQUIRE_OPENINVENTOR])
    AC_REQUIRE([mni_REQUIRE_QTGL])

    AC_LANG_PUSH(C++)
    mni_REQUIRE_LIB(SoQt,[#include <Inventor/Qt/SoQt.h>],[SoQt::mainLoop();])
    AC_LANG_POP
])

