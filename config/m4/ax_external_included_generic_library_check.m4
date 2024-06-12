# ===========================================================================
#  https://https://gite.lirmm.fr/doccy/doccy-dev-tools/autoconf
# ===========================================================================
#
# Serial 4
#
# SYNOPSIS
#
#   This macro defines options for using either system available or included
#   external library and computes the compiler and linker flags accordingly
#
#   AX_EXTERNAL_INCLUDED_GENERIC_LIBRARY_CHECK(C_CXX_FLAGS,
#                                              LIB_FLAGS,
#                                              NAME,
#                                              INCLUDED_LIB_CONFIG,
#                                              [LIBRARY_NAME],
#                                              [VERSION],
#                                              [FUNCTION],
#                                              [OTHER_LIBRARIES],
#                                              [MESSAGE-IF-NOT-FOUND])
#
# DESCRIPTION
#
#   The following configure script options are defined:
#   - '--with-included-NAME[=yes|no|check]' (default to is 'check')
#   - '--without-included-NAME' (alias for '--without-included-NAME=no').
#   - '--with-NAME-prefix=DIR'
#
#   TL;DR:
#   - The NAME_CPPFLAGS and NAME_LDFLAGS variables are defined with the
#     correct values to pass respectively to the compiler and to the linker.
#   - The INCLUDED_NAME automake conditional is defined to true iff the
#     included library is used.
#
#   This macro checks for library LIBRARY_NAME (if LIBRARY_NAME is empty,
#   then use NAME instead) according to the options '--with-included-NAME'
#   and '--with-NAME-prefix'.
#
#   If both options --with-NAME-prefix and --with-included-NAME[=yes] are
#   provided, configuration stops with an explicit error message.
#
#   If --with-NAME-prefix is set, uses the specified prefix to look for the
#   library.
#
#   If --with-included-NAME is not given or set to 'check' or 'no', then an
#   existing library is search on the system:
#   - If not found and --with-included-NAME is set to 'no', the configuration
#     stops with the MESSAGE-IF-NOT-FOUND (or a default one if empty).
#   - If not found and --with-included-NAME is set to 'check', then it acts
#     like if the --with-included-NAME was set to 'yes'.
#
#   If --with-included-NAME is set to 'yes' (default if argument value is not
#   explicitely set) or if the library hasn't been found on the system, then
#   the included library is used. In such case, The INCLUDED_LIB_CONFIG is a
#   comma separated list having between one and three values (if it contains
#   only one value, its is not required to quote the value but otherwise, the
#   list parameter must be quoted):
#   - The first (mandatory) element of this parameter is the relative path to
#     the included library,
#   - The second (optional) element is the C/CXX compiler flag(s) to set. If
#     emtpy, then append the `include` subdirectory to the included library
#     relative path (*i.e.*, `-I<included_lib_path>/include`),
#   - The third element is the linker flag(s) to set. If emtpy, then append
#     the `lib` subdirectory to the included library relative path and add
#     the library name to use (*i.e.*, `-L<included_lib_path>/lib
#     -l<libname>`, where <name> is computed using LIBRARY_NAME or NAME by
#     removing the `lib` prefix if any).
#
#   To detect system installed libraries, if the PKG_CONFIG shell variable is
#   not empty, it is assumed that it has been set to a valid path to program
#   `pkg-config` by a preceeding call to the macro PKG_PROG_PKG_CONFIG() and
#   then that the PKG_CHECK_MODULES() macro is available on the host system.
#   In such case, the library existance is first checked by using
#   PKG_CHECK_MODULES() (and the given VERSION constraint -- if any -- must
#   be satisfied). If the PKG-CONFIG shell variable is empty or if the
#   library is not found using this tool, then the FUNCTION is searched by
#   using the legacy AC_SEARCH_LIBS() macro (trying the LIBRARY_NAME first,
#   then the given OTHER_LIBRARIES if any).
#
#   If some buggy situation occurs, then a message is shown that proposes to
#   send an e-mail to the address specified in the PACKAGE_BUGREPORT shell
#   variable (or to the "package maintainers" if empty).
#
#   In the end, the given compiler and linker flags prefixed by NAME_ as well
#   as the INCLUDED_NAME automake conditional are set according to detection.
#
#   This macro is mainly devoted to the following macros:
#     - AX_EXTERNAL_INCLUDED_C_LIBRARY_DETECTION()
#     - AX_EXTERNAL_INCLUDED_CXX_LIBRARY_DETECTION()
#
#   ANY CHANGE IN THIS FILE MUST ENSURE NOT BREAKING THE ABOVE MACROS!
#
#   Note: This macro relies on the AX_REQUIRE_DEFINED(), the
#   AX_APPEND_COMPILE_FLAGS() and the AX_APPEND_LINK_FLAGS() external macros,
#   which are all available from the autoconf archive
#   (https://www.gnu.org/software/autoconf-archive/ax_require_defined.html,
#   https://www.gnu.org/software/autoconf-archive/ax_append_link_flags.html
#   and
#   https://www.gnu.org/software/autoconf-archive/ax_append_compile_flags.html).
#   It also relies to the AX_EXTERNAL_INCLUDED_LIBRARY_DECLARE_OPTIONS()
#   macro coming with this one.
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

AC_DEFUN([AX_EXTERNAL_INCLUDED_GENERIC_LIBRARY_CHECK],
[dnl
AX_REQUIRE_DEFINED([AX_APPEND_COMPILE_FLAGS])dnl
AX_REQUIRE_DEFINED([AX_APPEND_LINK_FLAGS])dnl
AX_REQUIRE_DEFINED([AX_EXTERNAL_INCLUDED_LIBRARY_DECLARE_OPTIONS])dnl

AX_EXTERNAL_INCLUDED_LIBRARY_DECLARE_OPTIONS([$3])
_priv_valid_found="no"

dnl Considers the library name is LIBRARY_NAME if not empty or NAME
dnl otherwise.
_priv_library=ifelse([$5], [], [$3], [$5])
dnl Computes the library name to link with the -l flag (the library name
dnl where the lib prefix is stripped if any).
_priv_library_name=${_priv_library#lib}

dnl In order to compute the included library parameters, we first define the
dnl given INCLUDED_LIB_CONFIG as a list.
m4_define([_INCLUDED_LIB_CONFIG], m4_dquote($4))
dnl The first element of INCLUDED_LIB_CONFIG list is the relative path.
_priv_included_path=m4_car(_INCLUDED_LIB_CONFIG)

dnl We thus remove the first element.
m4_define([_INCLUDED_LIB_CONFIG], m4_cdr(_INCLUDED_LIB_CONFIG))
dnl Now, the of the list (if not empty) is the C/C++ compiler flags.
_priv_default_c_cxx_flags="m4_car(_INCLUDED_LIB_CONFIG)"
dnl If the list is empty or if its head is the empty string, then use some
dnl default computed flags.
AS_IF([test -z "${_priv_default_c_cxx_flags}"],
       [_priv_default_c_cxx_flags="-I${ac_abs_confdir}/${_priv_included_path}/include"])dnl

dnl We remove the head of the list if not empty.
m4_define([_INCLUDED_LIB_CONFIG], m4_cdr(_INCLUDED_LIB_CONFIG))
dnl Now, the head of the list (if not empty) is the C/C++ linker flags.
_priv_default_ld_flags="m4_car(_INCLUDED_LIB_CONFIG)"
dnl If the list is empty or if its head is the empty string, then use some
dnl default computed flags.
AS_IF([test -z "${_priv_default_ld_flags}"],
      [_priv_default_ld_flags="-L${ac_pwd}${ac_dir_suffix}/${_priv_included_path}/lib -l${_priv_library_name}"])dnl

dnl We define the message to display if the pkg-config fails to find the
dnl library.
m4_define([_PKG_CONFIG_LIB_NOT_FOUND_MSG],
          [Library $3 ifelse([$6], [], [], [($6)]) not found using ${PKG_CONFIG} (with PKG_CONFIG_PATH='${PKG_CONFIG_PATH}'), falling back to the legacy library search.])dnl

dnl Check if some specific path is provided for NAME or not.
AS_IF([test "x${with_$3[]_prefix}" = "x"],
      [dnl Here, no specific path is provided to look for custom NAME
       dnl library.
       AS_IF([test "x${with_included_$3}" != "xyes"],
             [dnl Try to find an already existing version of NAME on system
              dnl with pkg-config.
              AS_IF([test "x${PKG_CONFIG}" != "x"],
                    [PKG_CHECK_MODULES([$3],
                                       [${_priv_library}],
                                       [_priv_CPPFLAGS=${$3[]_CFLAGS}
                                        _priv_LDFLAGS=${$3[]_LIBS}
                                        _priv_valid_found="yes"],
                                        [AC_MSG_NOTICE([_PKG_CONFIG_LIB_NOT_FOUND_MSG])])])
              AS_IF([test "x${_priv_valid_found}" = "xno"],
                    [dnl No valid system installed library found using
                     dnl pkg-config.
                     _priv_orig_LIBS="${LIBS}"
                     AC_SEARCH_LIBS([$7], [${_priv_library} $8],
                                    [_priv_valid_found="yes"])
                     _priv_LDFLAGS="${ac_cv_search_$7}"
                     LIBS="${_priv_orig_LIBS}"])
              AS_IF([test "x${_priv_valid_found}" = "xno"],
                    [dnl No valid system installed library found. Using
                     dnl legacy search.
                     AS_IF([test "x${with_included_$3}" = "xcheck"],
                           [dnl Fallback to the included library.
                            with_included_$3="yes"
                            AC_MSG_NOTICE([Library $3 not found on system, using the included one (${_priv_included_path}).])],
                           [dnl Aborting with an explicit error message
                            dnl (either the MESSAGE-IF-NOT-FOUND or the
                            dnl default error message.
                            AC_MSG_FAILURE([ifelse([$9], [],
                                           [

Unable to find required $3 library.

You can re-run the configure script by allowing the use of the
included $3 library.
],
                                           [$9])])])])])],
      [dnl Here some specific path for NAME is given. First ensure that there
       dnl is no conflicting use of --with-included-NAME option.
       AS_IF([test "x${with_included_$3}" = "xyes"],
             [AC_MSG_ERROR([

You can't use both --with-included-$3 and --with-$3-prefix options
])])
       dnl Check the library validity with pkg-config using the given path.
       AS_IF([test "x${PKG_CONFIG}" != "x"],
             [PKG_CONFIG_PATH="${with_$3[]_prefix}/lib/pkgconfig:${PKG_CONFIG_PATH}"
              export PKG_CONFIG_PATH
              PKG_CHECK_MODULES([$3],
                                [$6],
                                [_priv_CPPFLAGS=${$3[]_CFLAGS}
                                 _priv_LDFLAGS=${$3[]_LIBS}
                                 _priv_valid_found="yes"],
                                [AC_MSG_NOTICE([_PKG_CONFIG_LIB_NOT_FOUND_MSG])])])
       AS_IF([test "x${_priv_valid_found}" = "xno"],
             [dnl No valid system custom installed library found. Using
              dnl legacy search.
              _priv_CPPFLAGS="-I${with_$3[]_prefix}/include"
              _priv_LDFLAGS="-L${with_$3[]_prefix}/lib"
              _priv_orig_CPPFLAGS=${$1}
              _priv_orig_LDFLAGS=${$2}
              _priv_orig_LIBS="${LIBS}"
              AX_APPEND_COMPILE_FLAGS([${_priv_CPPFLAGS}], [$1])
              AX_APPEND_LINK_FLAGS([${_priv_LDFLAGS}], [$2])
              AC_SEARCH_LIBS([$7], [${_priv_library} $8],
                             [_priv_valid_found="yes"],
                             [AC_MSG_FAILURE([

No valid $3 library found using the path specified by option
--with-$3-prefix='${with_$3[]_prefix}'.
])])
              $1=${_priv_orig_CPPFLAGS}
              $2=${_priv_orig_LDFLAGS}
              _priv_LDFLAGS+="${ac_cv_search_$7}"
              LIBS="${_priv_orig_LIBS}"])])dnl

AS_IF([test "x${_priv_valid_found}" = "xno"],
      [dnl No library found on the system (either because option
       dnl --with-included-NAME[=yes] has been given or because there is no
       dnl valid installed library and --with-included-NAME was not set (or
       dnl was given the 'check' parameter).
      AS_IF([test "x${with_included_$3}" = "xyes"],
             [dnl Setting the flags for using the included library.
              _priv_CPPFLAGS="${_priv_default_c_cxx_flags}"
              _priv_LDFLAGS="${_priv_default_ld_flags}"
              m4_pushdef([AS_LITERAL_IF], [])# Undefine literal checking which issues a warning.
              AC_CONFIG_SUBDIRS([${_priv_included_path}])
              m4_popdef([AS_LITERAL_IF])# Restore the macro.
              _priv_valid_found="yes"],
             [dnl This case should not occur.
              _priv_valid_found="bug"])],
      [dnl Some valid installed library was found.
       AS_IF([test "x${with_included_$3}" = "xyes"],
             [dnl This case should not occur.
              _priv_valid_found="bug"],
             [dnl Update the value of the with_included_NAME variable to 'no'.
              with_included_$3="no"])])dnl

dnl If some bug occurs, we display an explicit message then abort
dnl configuration.
AS_IF([test "x${_priv_valid_found}" = "xbug"],
      [AC_MSG_FAILURE([

Something wrong occurs with the configure script and the $3 library detection.

Please send a bug report to ${PACKAGE_BUGREPORT:-package maintainers} describing at least the options passed to this configure script.
])])dnl

dnl Set the computed flags variables.
$3[]_$1=${_priv_CPPFLAGS}
$3[]_$2=${_priv_LDFLAGS}

dnl Set INCLUDED_NAME conditional for automake;
AM_CONDITIONAL([INCLUDED_$3], [test "x${with_included_$3}" = "xyes"])dnl

])dnl AX_EXTERNAL_INCLUDED_GENERIC_LIBRARY_CHECK
dnl Local variables:
dnl mode: autoconf
dnl fill-column: 77
dnl End:
