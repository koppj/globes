dnl Here follow some house-made macros to solve the SUSE libf2c problem

dnl Check for g2c
AC_DEFUN([AC_CHECK_G2C]),[
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
		F2CINC="\$(top_builddir)/libf2c"
	fi

	if test -n "$F2CCONVENIENCE"; then
		echo "libf2c/libg2c seems to be not installed in your system"
		echo "Using convenience library instead ..."
		f2c_convenience="yes"
	fi
else
	F2CCONVENIENCE="\$(top_builddir)/libf2c/libf2c.la"
	F2CINC="\$(top_builddir)/libf2c"
fi
dnl This control wether the stuff in `libf2c' is built
AM_CONDITIONAL(WANT_LIBF2C,test x$f2c_convenience = xyes)
dnl These defines are needed for proper compilation of the convenienece library
AC_DEFINE(USE_STRLEN)
AC_DEFINE(NON_ANSI_RW_MODES)
AC_DEFINE(UIOLEN_int)
AC_DEFINE(NON_UNIX_STDIO)
dnl This contains the correct linker flag for using the convience library
AC_SUBST(F2CCONVENIENCE)
AC_SUBST(F2CINC)
])

dnl Experimental support for making rpm's
dnl This is largely adapted from the rpm.m4 macro by
dnl Dale K. Hawkins <dhawkins@cdrgts.com>

dnl AM_RPM_INIT
dnl Figure out how to create rpms for this system and setup for an
dnl automake target

AC_DEFUN([AM_RPM_INIT],
[dnl
AC_REQUIRE([AC_CANONICAL_HOST])
dnl Find the RPM (rpmbuild) program
AC_ARG_WITH(rpm-prog,[  --with-rpm-prog=PROG   Which rpm to use (optional)],
            rpm_prog="$withval", rpm_prog="rpmbuild")

AC_ARG_ENABLE(rpm-rules, [  --enable-rpm-rules       Try to create rpm make rules [[default=no]]],
                enable_rpm_rules="$withval",enable_rpm_rules=no)

AC_ARG_WITH(rpm-extra-args, [  --with-rpm-extra-args=ARGS       Run rpm with extra arguments (defaults to none)],
                rpm_extra_args="$withval", rpm_extra_args="")

AC_ARG_WITH(rpm-targets, [  --with-rpm-targets=ARGS       Run rpm with targets (default is host_cpu)],
                rpm_targets="$withval", rpm_targets="$host_cpu")


dnl echo enable_rpm_rules is $enable_rpm_rules
dnl echo rpm_prog is $rpm_prog

RPM_ARCHS=$rpm_targets
AC_SUBST(RPM_ARCHS)

  RPM_TARGET=""

  if test x$enable_rpm_rules = xno ; then
     echo "Not trying to build rpms for your system (use --enable-rpm-rules to override) "
     no_rpm=yes
  else
    if test x$rpm_prog != x ; then
       if test x${RPM_PROG+set} != xset ; then
          RPM_PROG=$rpm_prog
       fi
    fi
dnl update to rpmbuild -- PH
dnl echo "RPM_PROG is $RPM_PROG"
dnl replace rpm with $rpm_prog in the next line -- PH
    AC_PATH_PROG(RPM_PROG, $rpm_prog, no)
    no_rpm=no
    if test "$RPM_PROG" = "no" ; then
echo "*** RPM Configuration Failed"
echo "*** Failed to find the rpmbuild program." 
echo "*** If you want to build rpm packages indicate the path to the rpmbuild"
echo "*** program using --with-rpm-prog=PROG"
echo "*** On older systems rpmbuild may be replaced by rpm itself."
echo "*** You can test this by calling rpm -ba, if rpm recognizes this option"
echo "*** rpm can be used instead of rpmbuild, just re-run configure using"
echo "*** --with-rpm-prog=rpm"
      no_rpm=yes
      RPM_MAKE_RULES=""
    else
      AC_MSG_CHECKING(how rpm sets %{_rpmdir})
      rpmdir=`rpm --eval %{_rpmdir}`
      if test x$rpmdir = x"%{_rpmdir}" ; then
        AC_MSG_RESULT([not set (cannot build rpms?)])
        echo *** Could not determine the value of %{_rpmdir}
        echo *** This could be because it is not set, or your version of rpm does not set it
        echo *** It must be set in order to generate the correct rpm generation commands
        echo ***
        echo *** You might still be able to create rpms, but I could not automate it for you
        echo *** BTW, if you know this is wrong, please help to improve the rpm.m4 module
        echo *** Send corrections, updates and fixes to dhawkins@cdrgts.com.  Thanks.
      else
        AC_MSG_RESULT([$rpmdir])
      fi
      AC_MSG_CHECKING(how rpm sets %{_rpmfilename})
	rpmfilename_list=""
	for arch in `echo $rpm_targets`
	do
      	rpmfilename_list="$rpmfilename_list $rpmdir/`rpm --eval %{_rpmfilename} | sed "s/%{ARCH}/${arch}/g" | sed "s/%{NAME}/$PACKAGE/g" | sed "s/%{VERSION}/${VERSION}/g" | sed "s/%{RELEASE}/${RPM_RELEASE}/g"`"
	done
      AC_MSG_RESULT([$rpmfilename_list])
AC_MSG_CHECKING(handling different build archs)	

dnl checking the rpm source directory where the tar-ball will end up...
     	AC_MSG_CHECKING(how rpm sets %{_sourcedir})
	rpmsourcedir=`rpm --eval %{_sourcedir}`
	AC_MSG_RESULT([$rpmsourcedir])	

	RPM_SOURCE_DIR=${rpmsourcedir}
	dnl checking the rpm buil directory where the tar-ball will end up...
     	AC_MSG_CHECKING(how rpm sets %{_builddir})
	rpmbuilddir=`rpm --eval %{_builddir}`
	AC_MSG_RESULT([$rpmbuilddir])	

	RPM_BUILD_DIR=${rpmbuilddir}
      RPM_DIR=${rpmdir}
      RPM_TARGET=$rpmfilename_list
      RPM_ARGS="-ba $rpm_extra_args"
      RPM_TARBALL=${PACKAGE}-${VERSION}.tar.gz
      RPM_VERSION=${PACKAGE}-${VERSION}
    fi
  fi

  case "${no_rpm}" in
    yes) make_rpms=false;;
    no) make_rpms=true;;
    *) AC_MSG_WARN([bad value ${no_rpm} for no_rpm (not making rpms)])
       make_rpms=false;;
  esac

AC_SUBST(RPM_SOURCE_DIR)
AC_SUBST(RPM_BUILD_DIR)	
  AC_SUBST(RPM_DIR)
  AC_SUBST(RPM_TARGET)
  AC_SUBST(RPM_ARGS)
  AC_SUBST(RPM_TARBALL)
 AC_SUBST(RPM_VERSION)


  RPM_CONFIGURE_ARGS=${ac_configure_args}
  AC_SUBST(RPM_CONFIGURE_ARGS)
])

