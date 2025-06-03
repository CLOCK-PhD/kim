dnl Process this file with autoconf to produce a configure script.

dnl Copyright © 2023-2025 -- IGH / LIRMM / CNRS / UM
dnl                          (Institut de Génétique Humaine /
dnl                          Laboratoire d'Informatique, de Robotique et de
dnl                          Microélectronique de Montpellier /
dnl                          Centre National de la Recherche Scientifique /
dnl                          Université de Montpellier)
dnl
dnl
dnl Auteurs/Authors:
dnl   - Rémy COSTA       <remy.costa@igh.cnrs.fr>
dnl   - William RITCHIE  <william.ritchie@igh.cnrs.fr>
dnl   - Alban MANCHERON  <alban.mancheron@lirmm.fr>
dnl
dnl
dnl Programmeurs/Programmers:
dnl   - Rémy COSTA       <remy.costa@igh.cnrs.fr>
dnl   - Alban MANCHERON  <alban.mancheron@lirmm.fr>
dnl
dnl
dnl Contact:
dnl   - KIM list         <kim@lirmm.fr>
dnl
dnl -------------------------------------------------------------------------
dnl
dnl Ce logiciel  est un  programme informatique  permettant  d'identifier des
dnl variations génomiques à partir de données brutes de séquençage.
dnl
dnl Ce logiciel est régi par la  licence CeCILL  soumise au droit français et
dnl respectant les principes  de diffusion des logiciels libres.  Vous pouvez
dnl utiliser, modifier et/ou redistribuer ce programme sous les conditions de
dnl la licence CeCILL telle que diffusée par  le CEA,  le CNRS et l'INRIA sur
dnl le site "http://www.cecill.info".
dnl
dnl En contrepartie de l'accessibilité au code source et des droits de copie,
dnl de modification et de redistribution accordés par cette licence, il n'est
dnl offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,
dnl seule une responsabilité  restreinte pèse  sur l'auteur du programme,  le
dnl titulaire des droits patrimoniaux et les concédants successifs.
dnl
dnl À  cet égard  l'attention de  l'utilisateur est  attirée sur  les risques
dnl associés  au chargement,  à  l'utilisation,  à  la modification  et/ou au
dnl développement  et à la reproduction du  logiciel par  l'utilisateur étant
dnl donné  sa spécificité  de logiciel libre,  qui peut le rendre  complexe à
dnl manipuler et qui le réserve donc à des développeurs et des professionnels
dnl avertis  possédant  des  connaissances  informatiques  approfondies.  Les
dnl utilisateurs  sont donc  invités  à  charger  et  tester  l'adéquation du
dnl logiciel  à leurs besoins  dans des conditions  permettant  d'assurer  la
dnl sécurité de leurs systêmes et ou de leurs données et,  plus généralement,
dnl à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.
dnl
dnl Le fait que  vous puissiez accéder  à cet en-tête signifie  que vous avez
dnl pris connaissance  de la licence CeCILL,  et que vous en avez accepté les
dnl termes.
dnl
dnl -------------------------------------------------------------------------
dnl
dnl This software is a computer program whose purpose is to indentify genomic
dnl from raw sequencing data.
dnl
dnl This software is governed by the CeCILL license under French law and
dnl abiding by the rules of distribution of free software. You can use,
dnl modify and/ or redistribute the software under the terms of the CeCILL
dnl license as circulated by CEA, CNRS and INRIA at the following URL
dnl "http://www.cecill.info".
dnl
dnl As a counterpart to the access to the source code and rights to copy,
dnl modify and redistribute granted by the license, users are provided only
dnl with a limited warranty and the software's author, the holder of the
dnl economic rights, and the successive licensors have only limited
dnl liability.
dnl
dnl In this respect, the user's attention is drawn to the risks associated
dnl with loading, using, modifying and/or developing or reproducing the
dnl software by the user in light of its specific status of free software,
dnl that may mean that it is complicated to manipulate, and that also
dnl therefore means that it is reserved for developers and experienced
dnl professionals having in-depth computer knowledge. Users are therefore
dnl encouraged to load and test the software's suitability as regards their
dnl requirements in conditions enabling the security of their systems and/or
dnl data to be ensured and, more generally, to use and operate it in the same
dnl conditions as regards security.
dnl
dnl The fact that you are presently reading this means that you have had
dnl knowledge of the CeCILL license and that you accept its terms.

AC_PREREQ([2.70])


dnl ###########################################################################
dnl #                           Project Informations                          #
dnl ###########################################################################


dnl ============================ Project Settings =============================

AX_GIT_VERSIONING_PROG_PATH([config])
AX_GIT_VERSIONING([config])

AC_INIT([kim],
        [GV_MAJOR.GV_MINOR.GV_MICRO],
        [kim@lirmm.fr],
        [kim],
        [https://github.com/CLOCK-PhD/kim])

AX_PACKAGE_INFOS(
  [PRETTY_NAME:k-mer Identification Metric.],
  [SHORT_DESCRIPTION:Kim is a k-mer identification metric program.],
  [FULL_DESCRIPTION:dnl
The kim program is an alignment-free software which allows to predict genomic dnl
variants from raw sequencing data.],
  [AUTHORS:Remy COSTA <remy.costa@igh.cnrs.fr>, William RITCHIE <william.ritchie@igh.cnrs.fr> and Alban MANCHERON  <alban.mancheron@lirmm.fr>],
  [LICENSE:CeCILL],
  [COPYRIGHT_DATE:2023-2025])dnl

AC_REVISION([GV_REVISION])

AC_DEFINE_UNQUOTED([BEGIN_KIM_NAMESPACE], [namespace kim {], [Starts the kim namespace.])
AC_DEFINE_UNQUOTED([END_KIM_NAMESPACE], [}], [Ends the kim namespace.])

AC_COPYRIGHT([
###############################################################################
#                                                                             #
#  Copyright © 2023-2025 -- IGH / LIRMM / CNRS / UM                           #
#                           (Institut de Génétique Humaine /                  #
#                           Laboratoire d'Informatique, de Robotique et de    #
#                           Microélectronique de Montpellier /                #
#                           Centre National de la Recherche Scientifique /    #
#                           Université de Montpellier)                        #
#                                                                             #
#                                                                             #
#  Auteurs/Authors:                                                           #
#    - Rémy COSTA       <remy.costa@igh.cnrs.fr>                              #
#    - William RITCHIE  <william.ritchie@igh.cnrs.fr>                         #
#    - Alban MANCHERON  <alban.mancheron@lirmm.fr>                            #
#                                                                             #
#                                                                             #
#  Programmeurs/Programmers:                                                  #
#    - Rémy COSTA       <remy.costa@igh.cnrs.fr>                              #
#    - Alban MANCHERON  <alban.mancheron@lirmm.fr>                            #
#                                                                             #
#                                                                             #
#  Contact:                                                                   #
#    - KIM list         <kim@lirmm.fr>                                        #
#                                                                             #
#  -------------------------------------------------------------------------  #
#                                                                             #
#  Ce logiciel  est un  programme informatique  permettant  d'identifier des  #
#  variations génomiques à partir de données brutes de séquençage.            #
#                                                                             #
#  Ce logiciel est régi par la  licence CeCILL  soumise au droit français et  #
#  respectant les principes  de diffusion des logiciels libres.  Vous pouvez  #
#  utiliser, modifier et/ou redistribuer ce programme sous les conditions de  #
#  la licence CeCILL telle que diffusée par  le CEA,  le CNRS et l'INRIA sur  #
#  le site "http://www.cecill.info".                                          #
#                                                                             #
#  En contrepartie de l'accessibilité au code source et des droits de copie,  #
#  de modification et de redistribution accordés par cette licence, il n'est  #
#  offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,  #
#  seule une responsabilité  restreinte pèse  sur l'auteur du programme,  le  #
#  titulaire des droits patrimoniaux et les concédants successifs.            #
#                                                                             #
#  À  cet égard  l'attention de  l'utilisateur est  attirée sur  les risques  #
#  associés au   chargement, à   l'utilisation, à  la modification  et/ou au  #
#  développement et   à la reproduction du  logiciel par l'utilisateur étant  #
#  donné  sa spécificité  de logiciel libre,  qui peut le rendre  complexe à  #
#  manipuler et qui le réserve donc à des développeurs et des professionnels  #
#  avertis  possédant  des  connaissances  informatiques  approfondies.  Les  #
#  utilisateurs sont   donc invités   à charger   et tester  l'adéquation du  #
#  logiciel à   leurs besoins  dans des  conditions permettant  d'assurer la  #
#  sécurité de leurs systèmes et ou de leurs données et, plus  généralement,  #
#  à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.         #
#                                                                             #
#  Le fait que   vous puissiez accéder à cet  en-tête signifie que vous avez  #
#  pris connaissance de la licence CeCILL,   et que vous en avez accepté les  #
#  termes.                                                                    #
#                                                                             #
#  -------------------------------------------------------------------------  #
#                                                                             #
#  This software is a computer program whose purpose is to indentify genomic  #
#  from raw sequencing data.                                                  #
#                                                                             #
#  This software is governed by the CeCILL license under French law and       #
#  abiding by the rules of distribution of free software. You can use,        #
#  modify and/ or redistribute the software under the terms of the CeCILL     #
#  license as circulated by CEA, CNRS and INRIA at the following URL          #
#  "http://www.cecill.info".                                                  #
#                                                                             #
#  As a counterpart to the access to the source code and rights to copy,      #
#  modify and redistribute granted by the license, users are provided only    #
#  with a limited warranty and the software's author, the holder of the       #
#  economic rights, and the successive licensors have only limited            #
#  liability.                                                                 #
#                                                                             #
#  In this respect, the user's attention is drawn to the risks associated     #
#  with loading, using, modifying and/or developing or reproducing the        #
#  software by the user in light of its specific status of free software,     #
#  that may mean that it is complicated to manipulate, and that also          #
#  therefore means that it is reserved for developers and experienced         #
#  professionals having in-depth computer knowledge. Users are therefore      #
#  encouraged to load and test the software's suitability as regards their    #
#  requirements in conditions enabling the security of their systems and/or   #
#  data to be ensured and, more generally, to use and operate it in the same  #
#  conditions as regards security.                                            #
#                                                                             #
#  The fact that you are presently reading this means that you have had       #
#  knowledge of the CeCILL license and that you accept its terms.             #
#                                                                             #
###############################################################################
])

AC_CONFIG_SRCDIR([src/kim.cpp])

AC_CONFIG_AUX_DIR([config])

dnl create a config.h file (Automake will add -DHAVE_CONFIG_H)
AC_CONFIG_HEADERS([config/config.h])
AH_TOP([
/******************************************************************************
*                                                                             *
*  Copyright © 2023-2025 -- IGH / LIRMM / CNRS / UM                           *
*                           (Institut de Génétique Humaine /                  *
*                           Laboratoire d'Informatique, de Robotique et de    *
*                           Microélectronique de Montpellier /                *
*                           Centre National de la Recherche Scientifique /    *
*                           Université de Montpellier)                        *
*                                                                             *
*                                                                             *
*  Auteurs/Authors:                                                           *
*    - Rémy COSTA       <remy.costa@igh.cnrs.fr>                              *
*    - William RITCHIE  <william.ritchie@igh.cnrs.fr>                         *
*    - Alban MANCHERON  <alban.mancheron@lirmm.fr>                            *
*                                                                             *
*                                                                             *
*  Programmeurs/Programmers:                                                  *
*    - Rémy COSTA       <remy.costa@igh.cnrs.fr>                              *
*    - Alban MANCHERON  <alban.mancheron@lirmm.fr>                            *
*                                                                             *
*                                                                             *
*  Contact:                                                                   *
*    - KIM list         <kim@lirmm.fr>                                        *
*                                                                             *
*  -------------------------------------------------------------------------  *
*                                                                             *
*  Ce logiciel  est un  programme informatique  permettant  d'identifier des  *
*  variations génomiques à partir de données brutes de séquençage.            *
*                                                                             *
*  Ce logiciel est régi par la  licence CeCILL  soumise au droit français et  *
*  respectant les principes  de diffusion des logiciels libres.  Vous pouvez  *
*  utiliser, modifier et/ou redistribuer ce programme sous les conditions de  *
*  la licence CeCILL telle que diffusée par  le CEA,  le CNRS et l'INRIA sur  *
*  le site "http://www.cecill.info".                                          *
*                                                                             *
*  En contrepartie de l'accessibilité au code source et des droits de copie,  *
*  de modification et de redistribution accordés par cette licence, il n'est  *
*  offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,  *
*  seule une responsabilité  restreinte pèse  sur l'auteur du programme,  le  *
*  titulaire des droits patrimoniaux et les concédants successifs.            *
*                                                                             *
*  À  cet égard  l'attention de  l'utilisateur est  attirée sur  les risques  *
*  associés au   chargement, à   l'utilisation, à  la modification  et/ou au  *
*  développement et   à la reproduction du  logiciel par l'utilisateur étant  *
*  donné  sa spécificité  de logiciel libre,  qui peut le rendre  complexe à  *
*  manipuler et qui le réserve donc à des développeurs et des professionnels  *
*  avertis  possédant  des  connaissances  informatiques  approfondies.  Les  *
*  utilisateurs sont   donc invités   à charger   et tester  l'adéquation du  *
*  logiciel à   leurs besoins  dans des  conditions permettant  d'assurer la  *
*  sécurité de leurs systèmes et ou de leurs données et, plus  généralement,  *
*  à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.         *
*                                                                             *
*  Le fait que   vous puissiez accéder à cet  en-tête signifie que vous avez  *
*  pris connaissance de la licence CeCILL,   et que vous en avez accepté les  *
*  termes.                                                                    *
*                                                                             *
*  -------------------------------------------------------------------------  *
*                                                                             *
*  This software is a computer program whose purpose is to indentify genomic  *
*  from raw sequencing data.                                                  *
*                                                                             *
*  This software is governed by the CeCILL license under French law and       *
*  abiding by the rules of distribution of free software. You can use,        *
*  modify and/ or redistribute the software under the terms of the CeCILL     *
*  license as circulated by CEA, CNRS and INRIA at the following URL          *
*  "http://www.cecill.info".                                                  *
*                                                                             *
*  As a counterpart to the access to the source code and rights to copy,      *
*  modify and redistribute granted by the license, users are provided only    *
*  with a limited warranty and the software's author, the holder of the       *
*  economic rights, and the successive licensors have only limited            *
*  liability.                                                                 *
*                                                                             *
*  In this respect, the user's attention is drawn to the risks associated     *
*  with loading, using, modifying and/or developing or reproducing the        *
*  software by the user in light of its specific status of free software,     *
*  that may mean that it is complicated to manipulate, and that also          *
*  therefore means that it is reserved for developers and experienced         *
*  professionals having in-depth computer knowledge. Users are therefore      *
*  encouraged to load and test the software's suitability as regards their    *
*  requirements in conditions enabling the security of their systems and/or   *
*  data to be ensured and, more generally, to use and operate it in the same  *
*  conditions as regards security.                                            *
*                                                                             *
*  The fact that you are presently reading this means that you have had       *
*  knowledge of the CeCILL license and that you accept its terms.             *
*                                                                             *
******************************************************************************/

#ifndef __KIM_CONFIG_H__
#define __KIM_CONFIG_H__
])
AH_BOTTOM([
#endif
])

AC_CONFIG_MACRO_DIR([config/m4])

dnl ========================= End of Project Settings =========================


dnl ###########################################################################
dnl #                Checking System and Configuration Options                #
dnl ###########################################################################


dnl ============================== System check ===============================

AC_CANONICAL_HOST

dnl Checks for programs.
AC_PROG_INSTALL
AC_PROG_CPP
AC_PROG_CC
AC_PROG_CXX
AC_PROG_MAKE_SET
AC_PROG_AWK
AC_PROG_LN_S
AC_PROG_EGREP
AC_PROG_SED
AC_PROG_LEX([noyywrap])
AX_PROG_BISON
BISON=${YACC}
AX_PROG_BISON_VERSION([3.6], [], [AC_MSG_WARN([
A (bison) parser generator version 3.6 (or greater) is needed if some built sources must be regenerated.
])
YACC=false
BISON=false
])
AS_IF([test "${LEX}" = ":"],
      [dnl
       LEX=false
       AC_MSG_WARN([
A ([f]lex) lexical analyser generator is needed if some built sources must be regenerated.
])])
PKG_PROG_PKG_CONFIG([0.25])

DEBIAN_PACKAGE_DESCRIPTION=`printf "${PACKAGE_FULL_DESCRIPTION}" | sed -e 's,^[[[:space:]]]*$,.,g' -e 's,^, ,g'`
AC_SUBST([DEBIAN_PACKAGE_DESCRIPTION])
AM_SUBST_NOTMAKE([DEBIAN_PACKAGE_DESCRIPTION])
AX_DEBIAN_PACKAGE([packaging/control],
                  [PRERM:packaging/prerm],
                  [POSTINST:packaging/postinst])

AC_LANG([C++])

dnl Ensure compiler is C++17 compatible
AX_CXX_COMPILE_STDCXX([17], [noext], [mandatory])

dnl Force some extra flags if not already present
AX_APPEND_COMPILE_FLAGS([-Wall -Wextra -pedantic], [CXXFLAGS])

dnl Using libtool
LT_INIT
LT_LANG([C++])
LT_PATH_LD
AX_COMPILER_FLAGS_LDFLAGS([EXTRA_LDFLAGS])
AX_APPEND_LINK_FLAGS([${EXTRA_LDFLAGS}], [LDFLAGS])

dnl Automake/Makefile settings
AM_INIT_AUTOMAKE([foreign no-texinfo.tex])
AM_SILENT_RULES([yes])

dnl =========================== End of System check ===========================


dnl ============================= Library Settings ============================

LIB_VERSION=GV_LIB_VERSION
LIB_VERSION_EXT=GV_LIB_VERSION_EXT

dnl LIB_VERSION is not the release number.
dnl Please read libtool documentation about versionning

LIBNAME=${PACKAGE_NAME#lib}
AC_SUBST([LIBNAME])
AC_SUBST([VERSION])
AC_SUBST([LIB_VERSION])
AC_SUBST([LIB_VERSION_EXT])
AC_DEFINE_UNQUOTED([LIB_VERSION],["GV_LIB_VERSION_EXT"],[Libtool version string])

dnl ========================= End of Library Settings =========================


dnl ============================= Code Statistics =============================

AX_CODE_STATISTICS()

dnl ========================= End of code statistics ==========================


dnl ================================= Doxygen =================================

AC_ARG_VAR([DX_DOXYGEN],
           [Documentation generator command (default: doxygen).])dnl

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
DX_INIT_DOXYGEN([${PACKAGE_PRETTY_NAME}], [doxygen.cfg], [docs])

AS_IF([test -z "${DX_DOXYGEN}"],
      [BUILD_DOC=false
       AC_MSG_WARN([
   The 'doxygen' program could not be found or documentation generation has
   been disable by providing the '--disable-doxygen-doc' option to the
   configure script.
   This isn't critical, but it means that you won't be able to create the
   documentation.
])],
      [AS_IF([test "${DX_FLAG_doc}" -ne 1],
             [BUILD_DOC=false
              AC_MSG_WARN([
   The documentation generation is disabled.
])],
             [BUILD_DOC=true
              DX_DOXYGEN_VERSION=`${DX_DOXYGEN} --version | cut -d' ' -f1`
              AX_COMPARE_VERSION([${DX_DOXYGEN_VERSION}],[lt],[${DX_MIN_VERSION}],
                                 [AC_MSG_WARN([
   Your version of the 'doxygen' program is quite old (v. ${DX_DOXYGEN_VERSION}).
   The generation of the documentation may fail or the result may look bad.
   This isn't critical, but it means that you won't be able to create a proper
   API documentation and thus to make correct releases.
])])])])

AM_CONDITIONAL([BUILD_DOC], [${BUILD_DOC}])
dnl Whatever the value of BUILD_DOC, the doc subdirectory need to be
dnl added to the configuration files see
dnl https://www.gnu.org/software/automake/manual/html_node/Unconfigured-Subdirectories.html
dnl for more details.
dnl The `doc/doxygen.cfg` and `doc/Makefile` targets are added to the `AC_CONFIG_FILES()` macro arguments at the end of this file
AM_CONDITIONAL([BUILD_DOC_MAN], [test ${DX_FLAG_man} = 1])

dnl ============================= End of Doxygen ==============================


dnl ================================== OpenMP =================================

dnl add option for compiling with OpenMP
AC_ARG_ENABLE([OpenMP],
              [AS_HELP_STRING([--disable-OpenMP],
                              [Disable the use of OpenMP])],
              [],
              [AX_OPENMP])

AX_APPEND_COMPILE_FLAGS([${OPENMP_CFLAGS}], [CFLAGS])
AX_APPEND_COMPILE_FLAGS([${OPENMP_CXXFLAGS}], [CXXFLAGS])
AX_APPEND_LINK_FLAGS([${OPENMP_LDFLAGS}], [LDFLAGS])

dnl ============================== End of OpenMP ==============================


dnl ============================== Code Coverage ==============================

dnl All formats of code coverage reports supported by the gcovr (until
dnl version 7.2) are enabled by default.

dnl Set the verbosity level of AX_GCOVR_CODE_COVERAGE() and associated
dnl macros to INFO (default is WARN).
GCOVR_AC_VERBOSITY(INFO)

dnl Apply a regex replacement in `*.gcov` files to merge files having
dnl different paths due to symlinks when built in some subdirectory.
dnl This is ugly, but it does the job...
GCOVR_GCOV_FILE_SED_COMMAND_ON_CUSTOM_BUILDDIR([s,0:Source:../\$(LIBNAME)/,0:Source:\$(top_srcdir)/../src/,])

dnl Disable TEXT/XML/JSON/SONARQUBE/COVERALLS/CSV/LEGACY code coverage report
dnl by default (by setting an empty target file). This will add
dnl `--enable-code-coverage-*format*[=target]` options in your
dnl configure script, for *format* in text, xml, sonarqube and
dnl coveralls.

dnl GCOVR_SUMMARY_TARGET()
GCOVR_TEXT_TARGET()
GCOVR_XML_TARGET()
dnl GCOVR_HTML_TARGET()
GCOVR_JSON_TARGET()
GCOVR_SONARQUBE_TARGET()
GCOVR_COVERALLS_TARGET()
GCOVR_CSV_TARGET()
dnl GCOVR_COBERTURA_TARGET()
GCOVR_JACOCO_TARGET()
GCOVR_CLOVER_TARGET()
GCOVR_LEGACY_HTML_TARGET()

dnl For the other *format*s the cnfigure script will be added the
dnl `--disable-code-coverage-*format*` option.

dnl Don't need to call the AX_CODE_COVERAGE() since it will be done by
dnl the AX_GCOVR_CODE_COVERAGE() macro
AX_GCOVR_CODE_COVERAGE()

dnl The macro AX_GCOVR_CODE_COVERAGE() will add the
dnl `--enable-code-coverage` and the `--with-gcov=GCOV` options (see
dnl AX_CODE_COVERAGE() macro documentation) and the new
dnl `--with-gcovr=GCOVR` option in configure script to provide a
dnl custom GCOVR program. The GCOVR environment variable can also be
dnl used to define the GCOVR program if not explicitly set.

dnl Finally, append the computed compiler/linker flags to enable the
dnl code coverage computation when code coverage is enabled.
AX_APPEND_COMPILE_FLAGS([${CODE_COVERAGE_CFLAGS}], [CFLAGS])
AX_APPEND_COMPILE_FLAGS([${CODE_COVERAGE_CPPFLAGS}], [CPPFLAGS])
AX_APPEND_COMPILE_FLAGS([${CODE_COVERAGE_CXXFLAGS}], [CXXFLAGS])
AX_APPEND_LINK_FLAGS([${CODE_COVERAGE_LIBS}], [LIBS])

dnl ========================== End of Code Coverage ===========================


dnl ================================== Debug ==================================

AX_C_CXX_ENABLE_DEBUG_OPTION([Passing argument 'print' does the same thing than 'yes' and
                              enables printing debugging informations on output too.
                              This option is implied by '--enable-code-coverage'],
                             [test "x${enable_code_coverage}" = "xyes"])

AS_CASE(["${enable_debug}"],
        [""], [enable_debug=no],
        [yes], [],
        [no],  [],
        [print], [],
        [AC_MSG_WARN([Value '${enable_debug}' is not valid for the '--enable-debug' option. Assuming 'yes'])
         enable_debug=yes])

dnl ============================== End of Debug ===============================


dnl ###########################################################################
dnl #                          Checking for libraries                         #
dnl ###########################################################################


dnl ================================= HTSlib ==================================

HTSlib_MIN_VERSION=1.17
AC_SUBST(HTSlib_MIN_VERSION)

AX_EXTERNAL_INCLUDED_CXX_LIBRARY_DETECTION([HTSlib],
                                           [[external/htslib],
                                            [-I${ac_pwd}/external/htslib/],
                                            [-L${ac_pwd}/external/htslib/ -lhts]],
                                           [htslib],
                                           [htslib > ${HTSlib_MIN_VERSION}],
                                           [bcf_strerror],
                                           [hts],
                                           [

Unable to find required HTSlib library (>= ${HTSlib_MIN_VERSION}).

You can re-run the configure script by allowing the use of the
included HTSlib library.

Alternatively, you can install this library using your package manager
(if up-to-date version is available):
# apt install libhts-dev # on Debian/Ubuntu-like distributions
# yum install htslib-devel # on Fedora/RedHat-like distributions

You can follow the instructions from
https://github.com/samtools/htslib/blob/master/INSTALL (see "System
Specific Details" section) if you prefer a custom installation and/or
a more recent version.

After installing this library, you have to run the configure script again.
])

dnl Don't use neither AX_APPEND_COMPILE_FLAGS() nor
dnl AX_APPEND_LINK_FLAGS() since:
dnl
dnl   1. The AX_EXTERNAL_INCLUDED_CPP_HEADER_DETECTION(), the
dnl      AX_EXTERNAL_INCLUDED_C_LIBRARY_DETECTION() and the
dnl      AX_EXTERNAL_INCLUDED_CXX_LIBRARY_DETECTION() macros already
dnl      check for validity when using already installed
dnl      library/headers.
dnl
dnl   2. If using the included library/headers, the
dnl      (pre)compiler/linker flags may not be correctly appended
dnl      since the library/header is most probably not already
dnl      available
dnl
dnl An acceptable option is to update the (pre)compiler/linker by
dnl appending the computed flags if no more test involving compilation
dnl is executed hereafter.
dnl
dnl A more secure (and recommended) strategy is to defer the
dnl compiler/linker flag updates at the end of the configure script,
dnl before generating the Makefiles using AC_CONFIG_COMMANDS_PRE()
AC_CONFIG_COMMANDS_PRE([CXXFLAGS="${CXXFLAGS} ${HTSlib_CXXFLAGS}"
                        LDFLAGS="${LDFLAGS} ${HTSlib_LDFLAGS}"])

dnl Disable releases if using the included HTSlib library.
AS_IF([test "x${with_included_HTSlib}" = "xyes"],
      [MAKE_RELEASE=false
       AC_MSG_WARN([
   You are using the included HTSlib library.
   This is fine for us, but it means that you won't be able to build releases.
])])

dnl ============================== End of HTSlib ==============================


dnl ###########################################################################
dnl #                         Checking for header files                       #
dnl ###########################################################################


dnl ================================== vcfpp ==================================

AX_EXTERNAL_INCLUDED_CPP_HEADER_DETECTION([vcfpp], [external/vcfpp], [vcfpp.h])

dnl Don't use neither AX_APPEND_COMPILE_FLAGS() nor
dnl AX_APPEND_LINK_FLAGS() since:
dnl
dnl   1. The AX_EXTERNAL_INCLUDED_CPP_HEADER_DETECTION(), the
dnl      AX_EXTERNAL_INCLUDED_C_LIBRARY_DETECTION() and the
dnl      AX_EXTERNAL_INCLUDED_CXX_LIBRARY_DETECTION() macros already
dnl      check for validity when using already installed
dnl      library/headers.
dnl
dnl   2. If using the included library/headers, the
dnl      (pre)compiler/linker flags may not be correctly appended
dnl      since the library/header is most probably not already
dnl      available
dnl
dnl An acceptable option is to update the (pre)compiler/linker by
dnl appending the computed flags if no more test involving compilation
dnl is executed hereafter.
dnl
dnl A more secure (and recommended) strategy is to defer the
dnl compiler/linker flag updates at the end of the configure script,
dnl before generating the Makefiles using AC_CONFIG_COMMANDS_PRE()
AC_CONFIG_COMMANDS_PRE([CPPFLAGS="${CPPFLAGS} ${vcfpp_CPPFLAGS}"])

dnl ============================== End of vcfpp ===============================


dnl ================================ BitMagic =================================

AX_EXTERNAL_INCLUDED_CPP_HEADER_DETECTION([BitMagic], [external/BitMagic/src], [bm.h])

dnl Don't use neither AX_APPEND_COMPILE_FLAGS() nor
dnl AX_APPEND_LINK_FLAGS() since:
dnl
dnl   1. The AX_EXTERNAL_INCLUDED_CPP_HEADER_DETECTION(), the
dnl      AX_EXTERNAL_INCLUDED_C_LIBRARY_DETECTION() and the
dnl      AX_EXTERNAL_INCLUDED_CXX_LIBRARY_DETECTION() macros already
dnl      check for validity when using already installed
dnl      library/headers.
dnl
dnl   2. If using the included library/headers, the
dnl      (pre)compiler/linker flags may not be correctly appended
dnl      since the library/header is most probably not already
dnl      available
dnl
dnl An acceptable option is to update the (pre)compiler/linker by
dnl appending the computed flags if no more test involving compilation
dnl is executed hereafter.
dnl
dnl A more secure (and recommended) strategy is to defer the
dnl compiler/linker flag updates at the end of the configure script,
dnl before generating the Makefiles using AC_CONFIG_COMMANDS_PRE()
AC_CONFIG_COMMANDS_PRE([CPPFLAGS="${CPPFLAGS} ${BitMagic_CPPFLAGS}"])

dnl ============================= End of BitMagic =============================


dnl ============================== optionparser ===============================

AX_EXTERNAL_INCLUDED_CPP_HEADER_DETECTION([optionparser], [external/optionparser], [optionparser.h])

dnl Don't use neither AX_APPEND_COMPILE_FLAGS() nor
dnl AX_APPEND_LINK_FLAGS() since:
dnl
dnl   1. The AX_EXTERNAL_INCLUDED_CPP_HEADER_DETECTION(), the
dnl      AX_EXTERNAL_INCLUDED_C_LIBRARY_DETECTION() and the
dnl      AX_EXTERNAL_INCLUDED_CXX_LIBRARY_DETECTION() macros already
dnl      check for validity when using already installed
dnl      library/headers.
dnl
dnl   2. If using the included library/headers, the
dnl      (pre)compiler/linker flags may not be correctly appended
dnl      since the library/header is most probably not already
dnl      available
dnl
dnl An acceptable option is to update the (pre)compiler/linker by
dnl appending the computed flags if no more test involving compilation
dnl is executed hereafter.
dnl
dnl A more secure (and recommended) strategy is to defer the
dnl compiler/linker flag updates at the end of the configure script,
dnl before generating the Makefiles using AC_CONFIG_COMMANDS_PRE()
AC_CONFIG_COMMANDS_PRE([CPPFLAGS="${CPPFLAGS} ${optionparser_CPPFLAGS}"])

dnl =========================== End of optionparser ===========================


dnl ================================ gzstream =================================

AC_SEARCH_LIBS([gzread], [z], [], [AC_MSG_ERROR([Library 'zlib' not found.])])
AC_SEARCH_LIBS([BZ2_bzRead], [bz2], [], [AC_MSG_ERROR([Library 'bz2' not found.])])


AX_EXTERNAL_INCLUDED_CPP_HEADER_DETECTION([gzstream], [external/gzstream/include], [gzstream.h])

dnl Don't use neither AX_APPEND_COMPILE_FLAGS() nor
dnl AX_APPEND_LINK_FLAGS() since:
dnl
dnl   1. The AX_EXTERNAL_INCLUDED_CPP_HEADER_DETECTION(), the
dnl      AX_EXTERNAL_INCLUDED_C_LIBRARY_DETECTION() and the
dnl      AX_EXTERNAL_INCLUDED_CXX_LIBRARY_DETECTION() macros already
dnl      check for validity when using already installed
dnl      library/headers.
dnl
dnl   2. If using the included library/headers, the
dnl      (pre)compiler/linker flags may not be correctly appended
dnl      since the library/header is most probably not already
dnl      available
dnl
dnl An acceptable option is to update the (pre)compiler/linker by
dnl appending the computed flags if no more test involving compilation
dnl is executed hereafter.
dnl
dnl A more secure (and recommended) strategy is to defer the
dnl compiler/linker flag updates at the end of the configure script,
dnl before generating the Makefiles using AC_CONFIG_COMMANDS_PRE()
AC_CONFIG_COMMANDS_PRE([CPPFLAGS="${CPPFLAGS} ${gzstream_CPPFLAGS}"])

dnl ============================= End of gzstream =============================


AC_CHECK_INCLUDES_DEFAULT

AC_HEADER_DIRENT
AC_HEADER_MAJOR
AC_FUNC_CLOSEDIR_VOID
AC_CHECK_HEADERS([\
    libgen.h \
    sys/stat.h sys/types.h \
    unistd.h \
  ], [], [AC_MSG_ERROR([Header file '${ac_header}' not found.])])

dnl Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T
AC_TYPE_UINT8_T
AC_TYPE_INT8_T
AC_TYPE_UINT16_T
AC_TYPE_INT16_T
AC_TYPE_UINT32_T
AC_TYPE_INT32_T
AC_TYPE_UINT64_T
AC_TYPE_INT64_T

dnl Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_CHECK_FUNCS([\
  memset mkdir rmdir \
  strerror strtol \
])


dnl ###########################################################################
dnl #               Producing files according to configuration                #
dnl ###########################################################################

AC_CONFIG_COMMANDS([fix-htslib],
                   [dnl
                    FIX_TAG="Autoconf Fix Makefile Targets when used as submodule"
                    cd ${srcdir}/external/htslib || continue
                    if test ! -f configure; then
                      if command -v autoreconf > /dev/null; then
                        AC_MSG_NOTICE([Generating HTSlib configuration script and install missing tools])
                        autoreconf -i
                      else
                        AC_MSG_ERROR([

In order to use the included HTSlib library, you need to install the 'autoreconf'
program on your system.
])
                      fi
                    fi
                    AC_MSG_NOTICE([Library HTSlib configuration script ready])
                    FILE_TOPATCH=Makefile
                    if test -f "${FILE_TOPATCH}"; then
                      if grep -q "${FIX_TAG}" "${FILE_TOPATCH}"; then
                        AC_MSG_NOTICE([File '${FILE_TOPATCH}' from HTSlib already patched])
                      else
                        sed -i 's,^distclean\( .*\)\?:,distclean-orig\1:,' "${FILE_TOPATCH}"
                        sed -i 's,^check\( .*\)\?:,check-disabled\1:,' "${FILE_TOPATCH}"
                        sed -i 's,\(-rm -f .*\)\bhtscodecs.mk\b,\1,' "${FILE_TOPATCH}"
                        sed -i 's,\btar -c\b.*|,( cd "\$(srcdir)" \&\& \$\${TAR-tar} chof - --exclude="autom4te.cache/*" --exclude=".git*" . ) |,' "${FILE_TOPATCH}"
                        sed -i 's,\btar -x\b,\$\${TAR-tar} xo,' "${FILE_TOPATCH}"
                        sed -i 's,^\(\t+cd \$(distdir) && \$(MAKE)\),\1 srcdir="\$(CURDIR)/\$(srcdir)",' "${FILE_TOPATCH}"
                        tab="$(echo $'\t')"
                        cat >> "${FILE_TOPATCH}" << FIX_TAG_EOF

### ${FIX_TAG} ###

MISSING_AUTOMAKE_TARGETS = dvi check installcheck

.PHONY: \$(MISSING_AUTOMAKE_TARGETS) distclean-orig check-disabled

\$(MISSING_AUTOMAKE_TARGETS):
$(echo $'\t')@echo "Fake target for '\$[]@' added by the kim configure script"

distclean: distclean-orig
${tab}test -n "\$(distdir)" \\
${tab}|| test '\$(srcdir)' = '.' \\
${tab}|| rm htslib_vars.mk htscodecs.mk htslib.mk Makefile

htslib_vars.mk:
${tab}@if test "\$(realpath \$(CURDIR))" = "\$(realpath \$(srcdir))"; then \\
${tab}  if test -e ".git"; then \\
${tab}    echo "Restoring '\$[]@' from git"; \\
${tab}    git restore "\$[]@"; \\
${tab}  else \\
${tab}    echo "Unable to restore '\$[]@'"; \\
${tab}    false; \\
${tab}  fi; \\
${tab}else \\
${tab}  echo "Rebuild '\$[]@' from source directory (\$(srcdir))"; \\
${tab}  cd "\$(srcdir)" && \$(MAKE) \$(AM_MAKEFLAGS) "\$[]@"; \\
${tab}  cp --force "\$(srcdir)/\$[]@" "\$[]@"; \\
${tab}fi

uninstall: uninstall-\$(SHLIB_FLAVOUR)
${tab}@echo "Target for '\$[]@' added by the kim configure script"
${tab}for f in \$(BUILT_PROGRAMS); do \\
${tab}  rm -f "\$(DESTDIR)\$(bindir)/\$\$[]f"; \\
${tab}done
${tab}if test -n "\$(BUILT_PLUGINS)"; then \\
${tab}  for f in \$(BUILT_PLUGINS); do \\
${tab}    rm -f "\$(DESTDIR)\$(plugindir)/\$\$[]f"; \\
${tab}  done; \\
${tab}fi
${tab}for f in "\$(DESTDIR)\$(includedir)/htslib"/*.h; do \\
${tab}  rm -f "\$\$[]f"; \\
${tab}done
${tab}rm -f "\$(DESTDIR)\$(libdir)/libhts.a"
${tab}for f in annot-tsv.1 bgzip.1 htsfile.1 tabix.1; do \\
${tab}  rm -f "\$(DESTDIR)\$(man1dir)/\$\$[]f"; \\
${tab}done
${tab}for f in faidx.5 sam.5 vcf.5; do \\
${tab}  rm -f "\$(DESTDIR)\$(man5dir)/\$\$[]f"; \\
${tab}done
${tab}rm -f "\$(DESTDIR)\$(man7dir)/htslib-s3-plugin.7"
${tab}rm -f "\$(DESTDIR)\$(pkgconfigdir)/htslib.pc"
${tab}rm -f external/htslib/htscodecs/htscodecs/version.h

uninstall-so:
${tab}@echo "Target for '\$[]@' added by the kim configure script"
${tab}rm -f "\$(DESTDIR)\$(libdir)/libhts.so.\$(PACKAGE_VERSION)"
${tab}rm -f "\$(DESTDIR)\$(libdir)/libhts.so"
${tab}rm -f "\$(DESTDIR)\$(libdir)/libhts.so.\$(LIBHTS_SOVERSION)"

uninstall-cygdll:
${tab}@echo "Target for '\$[]@' added by the kim configure script"
${tab}rm -f "\$(DESTDIR)\$(bindir)/cyghts-\$(LIBHTS_SOVERSION).dll"
${tab}rm -f "\$(DESTDIR)\$(libdir)/libhts.dll.a"

uninstall-dll:
${tab}@echo "Target for '\$[]@' added by the kim configure script"
${tab}rm -f "\$(DESTDIR)\$(bindir)/hts-\$(LIBHTS_SOVERSION).dll"
${tab}rm -f "\$(DESTDIR)\$(libdir)/hts.dll.a"

uninstall-dylib:
${tab}@echo "Target for '\$[]@' added by the kim configure script"
${tab}rm -f "\$(DESTDIR)\$(libdir)/libhts.\$(PACKAGE_VERSION).dylib"
${tab}rm -f "\$(DESTDIR)\$(libdir)/libhts.dylib"
${tab}rm -f "\$(DESTDIR)\$(libdir)/libhts.\$(LIBHTS_SOVERSION).dylib"

### End of ${FIX_TAG} ###
FIX_TAG_EOF
                      AC_MSG_NOTICE([File '${FILE_TOPATCH}' from HTSlib patched])
                      fi
                    else
                      AC_MSG_WARN([Unable to apply fix on HTSlib, some make targets may fail])
                    fi
                    cd - > /dev/null
                   ])

AC_CONFIG_FILES([
  Makefile
  coverage/Makefile
  doc/doxygen.cfg doc/mainpage.dox doc/Makefile
  packaging/kim.pc
  packaging/Makefile
  resources/Makefile
  src/Makefile
  tests/Makefile
])
AC_OUTPUT
