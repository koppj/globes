dnl Here follow some house-made macros to solve the SUSE libf2c problem

dnl Check for g2c
AC_DEFUN([AC_CHECK_G2C],[
AC_CHECK_LIB(g2c,main,[],[NOF2C="no"])
])

dnl Check for f2c and taking care of the SUSE libf2c problem 
dnl and alse taking into account that this thing is sometimes
dnl called g2c!
AC_DEFUN([AC_CHECK_F2C],[
NOF2C="yes"
dnl Check for libf2c and taking care of the SUSE libf2c problem 
AC_CHECK_LIB(f2c,main,[],[
globes_temp_lib="$LIBS"
LIBS="-lf2c $LIBS"
AC_CACHE_CHECK([for libf2c on SUSE],globes_cv_libf2c,
[AC_LINK_IFELSE([AC_LANG_PROGRAM([[int MAIN__;]], [[
int main (void) 
{ 
 exit(0);
}]])],[globes_cv_libf2c="yes"],[globes_cv_libf2c="no"])
])
if test $globes_cv_libf2c = yes; then
 AC_DEFINE(HAVE_LIBF2C)
else
 LIBS="$globes_temp_lib"
dnl nothing like libf2c found perhaps libg2c exists ...
 AC_CHECK_G2C
 AC_CHECK_HEADERS([g2c.h])
fi
])
])

dnl Check for blas and taking care of the SUSE libf2c problem 
AC_DEFUN([AC_CHECK_SUSE_blas],[
AC_CHECK_LIB(blas,main,[],[
globes_temp_lib="$LIBS"
LIBS="-lblas $LIBS"
AC_CACHE_CHECK([blas on SUSE],globes_cv_libblas,
[AC_LINK_IFELSE([AC_LANG_PROGRAM([[int MAIN__;]], [[
int main (void) 
{ 
 exit(0);
}]])],[globes_cv_libblas="yes"],[globes_cv_libblas="no"])
])
if test $globes_cv_libblas = yes; then
 AC_DEFINE(HAVE_LIBBLAS)
else
 LIBS=$globes_temp_lib && NOBLAS="no"
fi
])
])

dnl Check for lapack and taking care of the SUSE libf2c problem
AC_DEFUN([AC_CHECK_SUSE_lapack],[ 
AC_CHECK_LIB(lapack,main,[],[
globes_temp_lib="$LIBS"
LIBS="-llapack $LIBS"
AC_CACHE_CHECK([lapack on SUSE],globes_cv_liblapack,
[AC_LINK_IFELSE([AC_LANG_PROGRAM([[int MAIN__;]], [[
int main (void) 
{ 
 exit(0);
}]])],[globes_cv_liblapack="yes"],[globes_cv_liblapack="no"])
])
if test $globes_cv_liblapack = yes; then
 AC_DEFINE(HAVE_LIBLAPACK)	
else
 LIBS="$globes_temp_lib" && NOLAPACK="no"
fi
])
])

dnl FIXME caching would be really nice !
dnl This macro serves to make it possible to use a convenience version
dnl of the BLAS and LAPACK functions GLoBES needs (zgeev_). By
dnl default it looks for an installed BLAS/LAPACK and if found they
dnl will be used in linking the program. If not found the convenience
dnl version is built in the `lapack' subdirectory and this one is
dnl linked statically. If the option `--enable-lapack-convience' is given
dnl no search for a installed BLAS/LAPACK is performed and the convience
dnl library is used.
AC_DEFUN([AC_LAPACK_CONVENIENCE],[
AC_ARG_ENABLE(lapack-convenience,[  --enable-lapack-convenience   Using the BLAS and LAPACK convience functions],bl_convenience="yes",bl_convenience="no")
if test $bl_convenience = "no"; then	
	NOBLAS="yes"
	NOLAPACK="yes"
	AC_CHECK_SUSE_blas
	AC_CHECK_SUSE_lapack
	BLCONVENIENCE=""
	if test $NOBLAS = "no"; then
		BLCONVENIENCE="\$(top_builddir)/lapack/libglblapack.la"
	fi

	if test $NOLAPACK = "no"; then
		BLCONVENIENCE="\$(top_builddir)/lapack/libglblapack.la"
	fi

	if test -n "$BLCONVENIENCE"; then
		echo "LAPACK/BLAS seems to be not installed in your system"
		echo "Using convenience library instead ..."
		bl_convenience="yes"
	fi
else
	BLCONVENIENCE="\$(top_builddir)/lapack/libglblapack.la"
fi
dnl This control wether the stuff in `lapack' is built
AM_CONDITIONAL(WANT_LAPACK,test x$bl_convenience = xyes)
dnl This contains the correct linker flag for using the convience library
AC_SUBST(BLCONVENIENCE)
])

dnl This macro serves to make it possible to use a convenience version
dnl of libf2c. By
dnl default it looks for an installed libf2c (or libg2c) and if found they
dnl will be used in linking the program. If not found the convenience
dnl version is built in the `libf2c' subdirectory and this one is
dnl linked statically. If the option `--enable-lapack-convience' is given
dnl no search for a installed libf2c is performed and the convience
dnl library is used.
AC_DEFUN([AC_LIBF2C_CONVENIENCE],[
AC_ARG_ENABLE(libf2c-convenience,[  --enable-libf2c-convenience   Using the libf2c convenience library],f2c_convenience="yes",f2c_convenience="no")
if test $f2c_convenience = "no"; then	
	NOF2C="yes"
	AC_CHECK_F2C
	F2CCONVENIENCE=""
	F2CINC=""
	if test $NOF2C = "no"; then
		F2CCONVENIENCE="\$(top_builddir)/libf2c/libf2c.la"
		F2CINC="-I\$(top_builddir)/libf2c"
	fi

	if test -n "$F2CCONVENIENCE"; then
		echo "libf2c/libg2c seems to be not installed in your system"
		echo "Using convenience library instead ..."
		f2c_convenience="yes"
	fi
else
	F2CCONVENIENCE="\$(top_builddir)/libf2c/libf2c.la"
	F2CINC="-I\$(top_builddir)/libf2c"
fi
dnl This control wether the stuff in `libf2c' is built
AM_CONDITIONAL(WANT_LIBF2C,test x$f2c_convenience = xyes)
dnl These defines are needed for proper compilation of the convenienece library
AC_DEFINE(USE_STRLEN,[],[Needed by convenience libf2c])
AC_DEFINE(NON_ANSI_RW_MODES,[],[Needed by convenience libf2c])
AC_DEFINE(UIOLEN_int,[],[Needed by convenience libf2c])
AC_DEFINE(NON_UNIX_STDIO,[],[Needed by convenience libf2c])
dnl fixing the silly header installation policy on red hat
if test -z "$F2CINC" ; then
AC_CHECK_HEADER(f2c.h,[],[
F2CINC="-I\$(top_builddir)/libf2c"
])
fi

dnl This contains the correct linker flag for using the convience library
AC_SUBST(F2CCONVENIENCE)
AC_SUBST(F2CINC)
])