dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.64)

AC_INIT([kim],
        [0.0.0],
        [remy.costa@igh.cnrs.fr],
        [kim],
        [https://github.com/CLOCK-PhD/kim])

PACKAGE_DESCRIPTION="Kim is a k-mer identification metric"
AC_SUBST([PACKAGE_DESCRIPTION])

AC_CONFIG_SRCDIR([src/kim.cpp])

AC_CONFIG_AUX_DIR([config])
dnl create a config.h file (Automake will add -DHAVE_CONFIG_H)
AM_CONFIG_HEADER([config/config.h])
AH_TOP([
#ifndef __KIM_CONFIG_H__
#define __KIM_CONFIG_H__
])
AH_BOTTOM([
#endif
])


dnl Every other copy of the package version number gets its value from here
AM_INIT_AUTOMAKE([foreign no-texinfo.tex])
AM_SILENT_RULES([yes])

AC_CONFIG_MACRO_DIR([config/m4])

AC_CANONICAL_HOST



dnl Checks for programs.
AC_PROG_INSTALL
AC_PROG_CC
AC_PROG_CXX

AC_PATH_PROG([CLOC], [cloc])
AC_PATH_PROG([OHCOUNT], [ohcount])

dnl ================================= Doxygen =================================

AC_ARG_VAR([DX_DOXYGEN], [Documenation generator command (default: doxygen).])
AC_PATH_PROG([DX_DOXYGEN], [${DX_DOXYGEN:-doxygen}])

DX_MIN_VERSION=1.9.6

DX_DOXYGEN_FEATURE(ON)
DX_DOT_FEATURE(ON)
DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(ON)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(ON)
DX_PS_FEATURE(OFF)
DX_INIT_DOXYGEN([${PACKAGE_NAME}], [doxygen.cfg], [docs])

AS_IF([test -z "${DX_DOXYGEN}"],
      [BUILD_DOC=false
       AC_MSG_WARN([The 'doxygen' program cound not be found. This isn't critical, but it means that you won't be able to create the documentation.])],
      [BUILD_DOC=true
      DX_DOXYGEN_VERSION=`${DX_DOXYGEN} --version | cut -d' ' -f1`
      AX_COMPARE_VERSION([${DX_DOXYGEN_VERSION}],[lt],[${DX_MIN_VERSION}],
                         [AC_MSG_WARN([Your version of the 'doxygen' program is quite old (v. ${DX_DOXYGEN_VERSION}) and the generation of the documentation may fail or the result mays look bad. This isn't critical, but it means that you won't be able to create a proper API documentation.])])])

AM_CONDITIONAL([BUILD_DOC], [${BUILD_DOC}])
AM_COND_IF([BUILD_DOC], [AC_CONFIG_FILES([doc/doxygen.cfg doc/Makefile])])
AM_CONDITIONAL([BUILD_DOC_MAN], [test ${DX_FLAG_man} = 1])

dnl ============================= End of Doxygen ==============================

ISODATE=`date +%Y-%m-%d`
AC_SUBST([ISODATE])

dnl Checks for libraries.

dnl Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([\
  unistd.h \
  sys/param.h sys/time.h \
  time.h \
  sys/mkdev.h sys/sysmacros.h \
  string.h memory.h fcntl.h dirent.h \
  sys/ndir.h \
  ndir.h alloca.h locale.h])

AC_HEADER_MAJOR
AC_FUNC_ALLOCA
AC_STRUCT_TM
AC_STRUCT_ST_BLOCKS
AC_FUNC_CLOSEDIR_VOID
AC_CHECK_FUNCS(mkfifo)
AC_CHECK_FUNC(mknod)

dnl Checks for typedefs, structures, and compiler characteristics.

dnl Checks for library functions.

AC_OUTPUT([
  Makefile src/Makefile
])

