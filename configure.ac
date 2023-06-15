dnl Process this file with autoconf to produce a configure script.

dnl Copyright © 2023      -- IGH / LIRMM / CNRS / UM
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

AC_PREREQ([2.64])


AC_INIT([kim],
        [0.0.0],
        [remy.costa@igh.cnrs.fr],
        [kim],
        [https://github.com/CLOCK-PhD/kim])

PACKAGE_DESCRIPTION="Kim is a k-mer identification metric"
AC_SUBST([PACKAGE_DESCRIPTION])


AC_COPYRIGHT([
###############################################################################
#                                                                             #
#  Copyright © 2023      -- IGH / LIRMM / CNRS / UM                           #
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
AM_CONFIG_HEADER([config/config.h])
AH_TOP([
/******************************************************************************
*                                                                             *
*  Copyright © 2023      -- IGH / LIRMM / CNRS / UM                           *
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

