# ===========================================================================
#  https://https://gite.lirmm.fr/doccy-dev-tools/autoconf
# ===========================================================================
#
# Serial 5
#
# SYNOPSIS
#
#   This macro defines options for using either system available or included
#   external headers and set the compiler flags accordingly.
#
#   AX_EXTERNAL_INCLUDED_CPP_HEADER_DETECTION(NAME,
#                                             INCLUDED_PATH,
#                                             FILE,
#                                             [MESSAGE-IF-NOT-FOUND])
#
# DESCRIPTION
#
#   The following configure script options are defined:
#   - `--with-included-NAME[=yes|no|check]` (default to is `check`)
#   - `--without-included-NAME` (alias for `--without-included-NAME=no`).
#   - `--with-NAME-path=DIR`
#
#   TL;DR:
#   - The `NAME_CPPFLAGS` variable is defined with the correct values to pass
#     to the compiler.
#   - The `INCLUDED_NAME` automake conditional is defined to true iff the
#     included headers are used.
#
#   This macro checks for library FILE according to the options
#   `--with-included-NAME` and `--with-NAME-path`.
#
#   If both options `--with-NAME-path` and `--with-included-NAME[=yes]` are
#   provided, configuration stops with an explicit error message.
#
#   If `--with-NAME-path` is set, uses the specified path to look for the
#   header.
#
#   If `--with-included-NAME` is not given or set to `check` or `no`, then an
#   existing header is search on the system:
#   - If not found and `--with-included-NAME` is set to `no`, the
#     configuration stops with the `MESSAGE-IF-NOT-FOUND` (or a default one
#     if empty).
#   - If not found and `--with-included-NAME` is set to `check`, then it acts
#     like if the `--with-included-NAME` was set to `yes`.
#
#   If `--with-included-NAME` is set to `yes` (default if argument value is
#   not explcitely set), then the included header path is used (expecting
#   that `INCLUDED_PATH/FILE` exists).
#
#   The header file detection on the host system is achieved by using the
#   `AC_CHECK_HEADERS()` macro.
#
#   If some buggy situation occurs, then a message is shown that proposes to
#   send an e-mail to the address specified in the `PACKAGE_BUGREPORT` shell
#   variable (or to the "package maintainers" if empty).
#
#   In the end, the `NAME_CPPFLAGS` compiler flag and the `INCLUDED_NAME`
#   automake conditional are set according to detection.
#
#   Notice that, accordingly to the well documented automake [Conditional
#   Subdirectories](https://www.gnu.org/software/automake/manual/automake.html#Conditional-Subdirectories)
#   section, the included header should be added to the subdirectories to
#   include in the distribution, even if the `--without-included-NAME` option
#   is given to the `configure` script. Although, you should add the included
#   header path to the automake `EXTRA_DIST` variable in your top level
#   `Makefile.am`:
#
#   ```automake
#   EXTRA_DIST = COPYING INSTALL NEWS
#   
#   EXTRA_DIST += INCLUDED_PATH
#   ```
#
#   Note: This macro relies on the `AX_REQUIRE_DEFINED()` and
#   `AX_APPEND_COMPILE_FLAGS()` external macros, which are both available
#   from the autoconf archive
#   (https://www.gnu.org/software/autoconf-archive/ax_require_defined.html
#   and
#   https://www.gnu.org/software/autoconf-archive/ax_append_compile_flags.html).
#   It also relies to the `AX_EXTERNAL_INCLUDED_LIBRARY_DECLARE_OPTIONS()`
#   macro coming with this one.
#
#   Remark: Someone might legitimately wonder why we prefer to use system
#   header files rather than those provided in the program. The motivation is
#   to avoid installing the included files if they are already present (in
#   the same version -- or worse -- in a different version) on the system and
#   thus to avoid conflicts and/or side effects for other programs.
#
# LICENSE
#
#   Copyright © 2015-2024 -- LIRMM / CNRS / UM
#                            (Laboratoire d'Informatique, de Robotique et de
#                            Microélectronique de Montpellier /
#                            Centre National de la Recherche Scientifique /
#                            Université de Montpellier)
#
#
#   Auteur/Author:
#     - Alban MANCHERON  <alban.mancheron@lirmm.fr>
#
#   Programmeurs/Programmers:
#     - Alban MANCHERON  <alban.mancheron@lirmm.fr>
#
#   -------------------------------------------------------------------------
#
#   Ce logiciel est régi par la  licence CeCILL  soumise au droit français et
#   respectant les principes  de diffusion des logiciels libres.  Vous pouvez
#   utiliser, modifier et/ou redistribuer ce programme sous les conditions de
#   la licence CeCILL telle que diffusée par  le CEA,  le CNRS et l'INRIA sur
#   le site "http://www.cecill.info".
#
#   En contrepartie de l'accessibilité au code source et des droits de copie,
#   de modification et de redistribution accordés par cette licence, il n'est
#   offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,
#   seule une responsabilité  restreinte pèse  sur l'auteur du programme,  le
#   titulaire des droits patrimoniaux et les concédants successifs.
#
#   À  cet égard  l'attention de  l'utilisateur est  attirée sur  les risques
#   associés  au chargement,  à  l'utilisation,  à  la modification  et/ou au
#   développement  et à la reproduction du  logiciel par  l'utilisateur étant
#   donné  sa spécificité  de logiciel libre,  qui peut le rendre  complexe à
#   manipuler et qui le réserve donc à des développeurs et des professionnels
#   avertis  possédant  des  connaissances  informatiques  approfondies.  Les
#   utilisateurs  sont donc  invités  à  charger  et  tester  l'adéquation du
#   logiciel  à leurs besoins  dans des conditions  permettant  d'assurer  la
#   sécurité de leurs systêmes et ou de leurs données et,  plus généralement,
#   à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.
#
#   Le fait que  vous puissiez accéder  à cet en-tête signifie  que vous avez
#   pris connaissance  de la licence CeCILL,  et que vous en avez accepté les
#   termes.
#
#   -------------------------------------------------------------------------
#
#   This software is governed by the CeCILL license under French law and
#   abiding by the rules of distribution of free software. You can use,
#   modify and/ or redistribute the software under the terms of the CeCILL
#   license as circulated by CEA, CNRS and INRIA at the following URL
#   "http://www.cecill.info".
#
#   As a counterpart to the access to the source code and rights to copy,
#   modify and redistribute granted by the license, users are provided only
#   with a limited warranty and the software's author, the holder of the
#   economic rights, and the successive licensors have only limited
#   liability.
#
#   In this respect, the user's attention is drawn to the risks associated
#   with loading, using, modifying and/or developing or reproducing the
#   software by the user in light of its specific status of free software,
#   that may mean that it is complicated to manipulate, and that also
#   therefore means that it is reserved for developers and experienced
#   professionals having in-depth computer knowledge. Users are therefore
#   encouraged to load and test the software's suitability as regards their
#   requirements in conditions enabling the security of their systems and/or
#   data to be ensured and, more generally, to use and operate it in the same
#   conditions as regards security.
#
#   The fact that you are presently reading this means that you have had
#   knowledge of the CeCILL license and that you accept its terms.

AC_DEFUN([AX_EXTERNAL_INCLUDED_CPP_HEADER_DETECTION],
[dnl
AC_REQUIRE_CPP()dnl
AX_REQUIRE_DEFINED([AX_APPEND_COMPILE_FLAGS])dnl
AX_REQUIRE_DEFINED([AX_EXTERNAL_INCLUDED_LIBRARY_DECLARE_OPTIONS])dnl

AX_EXTERNAL_INCLUDED_LIBRARY_DECLARE_OPTIONS([$1], [true])
valid_$1_found="no"

dnl Check if some specific path is provided for NAME or not.
AS_IF([test -z "${with_$1_path}"],
      [dnl Here, no specific path is provided to look for custom NAME header.
       AS_IF([test "${with_included_$1}" != "yes"],
             [dnl Try to find an already existing version of NAME on system.
              AC_CHECK_HEADERS([$3], [valid_$1_found="yes"])
              AS_IF([test "${valid_$1_found}" = "no"],
                    [dnl No valid system installed header found.
                     AS_IF([test "${with_included_$1}" = "check"],
                           [dnl Fallback to the included header.
                            with_included_$1="yes"
                            AC_MSG_NOTICE([Header(s) of $1 not found on system, using the one(s) included in '$2'.])],
                           [dnl Aborting with an explicit error message
                            dnl (either the MESSAGE-IF-NOT-FOUND or the
                            dnl default error message.
                            AC_MSG_FAILURE([m4_default([$4],
                                           [

Unable to find required $1 header(s).

You can re-run the configure script by allowing the use of the
included $1 header(s).
])])])])])],
      [dnl Here some specific path for NAME is given. First ensure that there
       dnl is no conflicting use of --with-included-NAME option.
       AS_IF([test "${with_included_$1}" = "yes"],
             [AC_MSG_ERROR([

You can't use both --with-included-$1 and --with-$1-path options
])])
       dnl Try to find an existing version of NAME on system using the given
       dnl path.
       CPPFLAGS_tmp="-I${with_$1_path}"
       CPPFLAGS_orig="${CPPFLAGS}"
       AX_APPEND_COMPILE_FLAGS([${CPPFLAGS_tmp}], [CPPFLAGS])
       AC_CHECK_HEADERS([$3],
                        [valid_$1_found="yes"],
                        [AC_MSG_FAILURE([

No valid $1 header found using the path specified by option
--with-$1-path='${with_$1_path}'.
])])
       CPPFLAGS=${CPPFLAGS_orig}])dnl

AS_IF([test "${valid_$1_found}" = "no"],
      [dnl No header found on the system (either because option
       dnl --with-included-NAME[=yes] has been given or because there is no
       dnl valid installed header and --with-included-NAME was not set (or
       dnl was given the 'check' parameter).
     AS_IF([test "${with_included_$1}" = "yes"],
             [dnl Setting the flags for using the included header.
              CPPFLAGS_tmp="-I${ac_abs_confdir}/$2"
              valid_$1_found="yes"],
             [dnl This case should not occur.
              valid_$1_found="bug"])],
      [dnl Some valid installed header was found.
       AS_IF([test "${with_included_$1}" = "yes"],
             [dnl This case should not occur.
              valid_$1_found="bug"],
             [dnl Update the value of the with_included_NAME variable to
              dnl 'no'.
              with_included_$1="no"])])dnl

dnl If some bug occurs, we display an explicit message then abort
dnl configuration.
AS_IF([test "${valid_$1_found}" = "bug"],
      [AC_MSG_FAILURE([

Something wrong occurs with the configure script and the $1 header detection.

Please send a bug report to ${PACKAGE_BUGREPORT:-package maintainers} describing at least the options passed to this configure script.
])])dnl

dnl Set the computed flags variable.
$1_CPPFLAGS=${CPPFLAGS_tmp}

dnl Set INCLUDED_NAME conditional for automake.
AM_CONDITIONAL([INCLUDED_$1], [test "${with_included_$1}" = "yes"])dnl

])dnl AX_EXTERNAL_INCLUDED_CPP_HEADER_DETECTION
dnl Local variables:
dnl mode: autoconf
dnl fill-column: 77
dnl End:
