dnl @synopsis mni_CXX_HAVE_KOENIG_LOOKUP
dnl
dnl Define CXX_HAVE_KOENIG_LOOKUP if the C++ compiler has 
dnl argument-dependent name lookup (a.k.a. Koenig lookup).
dnl
dnl @version $Id: mni_cxx_have_koenig_lookup.m4,v 1.2 2001/09/11 02:15:39 stever Exp $
dnl @author Steve Robbins
dnl
AC_DEFUN([mni_CXX_HAVE_KOENIG_LOOKUP],
    [AC_CACHE_CHECK(whether the compiler implements Koenig lookup,
                    ac_cv_cxx_have_koenig_lookup,
                    [AC_LANG_PUSH(C++)
                     AC_TRY_COMPILE([
    namespace N1 { 
    	class C {}; 
    	void f1(const C& c) {} 
    }

    namespace N2 {
    	void f2() {
	    N1::C x;
	    f1(x);     // resolves to N1::f1() if we have Koenig lookup,
                       // otherwise this will fail to compile.
        }
    }
    ],[],
		     ac_cv_cxx_have_koenig_lookup=yes,
                     ac_cv_cxx_have_koenig_lookup=no)
                     AC_LANG_POP])
    if test "$ac_cv_cxx_have_koenig_lookup" = yes; then
  	AC_DEFINE(CXX_HAVE_KOENIG_LOOKUP,1,
                  [define to 1 if the compiler implements Koenig lookup])
    fi
])
