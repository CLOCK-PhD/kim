###############################################################################
#                                                                             #
#  Copyright © 2023-2024 -- IGH / LIRMM / CNRS / UM                           #
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

AM_CXXFLAGS =				\
  -I@top_srcdir@			\
  -DCXX="\"$(CXX)\""			\
  -DARCH="\"$(build_cpu)\""		\
  -DOS="\"$(build_os)\""		\
  -DPACKAGE_DATADIR="\"$(pkgdatadir)\""

library_includedir = @includedir@/@PACKAGE@-@VERSION@
library_include_HEADERS =	\
  bounded_size_string.h		\
  fastq_file_reader.h		\
  file_reader.h			\
  kim.h				\
  kim_exception.h		\
  kim_settings.h		\
  kmer_nodes_subindex.h		\
  kmer_variant_edges_subindex.h	\
  kmer_variant_graph.h		\
  sort_helper.h                 \
  variant_identification.h	\
  variant_nodes_index.h

# Bug fix for path with spaces
bindir="@bindir@"
libdir="@libdir@"
# End of bug fix

lib_LTLIBRARIES = libkim.la
check_LTLIBRARIES = libkim-debug.la

libkim_ladir = @abs_srcdir@
libkim_la_SOURCES =						\
  bounded_size_string.h bounded_size_string.cpp			\
  dna_file_reader.h dna_file_reader.cpp				\
  fastq_file_reader.h fastq_file_reader.cpp			\
  file_reader.h file_reader.cpp					\
  kim_exception.h						\
  kim_settings.h kim_settings.cpp				\
  kmer_nodes_subindex.h kmer_nodes_subindex.cpp			\
  kmer_variant_edges_subindex.h kmer_variant_edges_subindex.cpp	\
  kmer_variant_graph.h kmer_variant_graph.cpp			\
  sort_helper.h sort_helper.cpp					\
  variant_identification.h variant_identification.cpp		\
  variant_nodes_index.h variant_nodes_index.cpp

libkim_debug_ladir = @abs_srcdir@
libkim_debug_la_SOURCES = $(libkim_la_SOURCES)
libkim_debug_la_CXXFLAGS = $(AM_CXXFLAGS) -UNDEBUG

GCDA_FILES = $(patsubst %.cpp,%.gcda,$(filter %.cpp,$(SOURCES)))
GCNO_FILES = $(GCDA_FILES:.gcda=.gcno)
GCOV_FILES = $(GCDA_FILES:.gcda=.gcov)

libkim_la_configdir   = @libdir@/@PACKAGE@-@VERSION@/include
libkim_la_config_DATA = $(CONFIG_HEADER)
libkim_la_LDFLAGS     = -version-info "@LIB_VERSION@" $(AM_LDFLAGS)

bin_PROGRAMS = kim
kim_SOURCES  = kim.h kim.cpp optionparser.h
kim_LDADD    = libkim.la

install-data-hook:
	$(AM_V_at)test "`readlink "$(DESTDIR)@includedir@/@PACKAGE@"`" == "@PACKAGE@-@VERSION@" \
	|| ( cd "$(DESTDIR)@includedir@" && $(LN_S) -f -n "@PACKAGE@-@VERSION@" "@PACKAGE@" )
	$(AM_V_at)test "`readlink "$(DESTDIR)@libdir@/@PACKAGE@"`" == "@PACKAGE@-@VERSION@" \
	|| ( cd "$(DESTDIR)@libdir@" && $(LN_S) -f -n "@PACKAGE@-@VERSION@" "@PACKAGE@" )

uninstall-hook:
	$(AM_V_at)test "`readlink "$(DESTDIR)@includedir@/@PACKAGE@"`" != "@PACKAGE@-@VERSION@" \
        || rm -rf "$(DESTDIR)@includedir@"/@PACKAGE@{,-@VERSION@}
	$(AM_V_at)test "`readlink "$(DESTDIR)@libdir@/@PACKAGE@"`" != "@PACKAGE@-@VERSION@" \
        || rm -rf "$(DESTDIR)@libdir@"/@PACKAGE@{,-@VERSION@}

MOSTLYCLEANFILES=*~
CLEANFILES=*~ $(GCDA_FILES) $(GCNO_FILES) $(GCOV_FILES)
DISTCLEANFILES=*~
MAINTAINERCLEANFILES=*~

