AC_INIT([RMINC], 0.4)

dnl Now find the compiler and compiler flags to use
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi

CC=`"${R_HOME}/bin/R" CMD config CC`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`

with_minc2=yes

smr_WITH_BUILD_PATH
mni_REQUIRE_MINC

dnl substitute CPPFLAGS and LIBS
AC_SUBST(CPPFLAGS)
AC_SUBST(LIBS)
dnl and do subsitution in the src/Makevars.in
AC_OUTPUT(src/Makevars)
