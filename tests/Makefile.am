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


########################
# Global configuration #
########################

AM_CXXFLAGS =				\
  -I$(top_srcdir)/src			\
  -I$(top_builddir)/src			\
  -DCXX="\"$(CXX)\""			\
  -DARCH="\"$(build_cpu)\""		\
  -DOS="\"$(build_os)\""		\
  -DPACKAGE_DATADIR="\"$(pkgdatadir)\""	\
  -DSRCDIR="\"$(abs_srcdir)\""

LDADD = ../src/libkim-debug.la

if INCLUDED_HTSlib
export LD_LIBRARY_PATH += :$(top_builddir)/external/htslib/
endif

check_PROGRAMS =
TESTS =


#############################################
# kim::BoundedSizeString class test program #
#############################################

check_PROGRAMS += test_bounded_size_string
TESTS          += test_bounded_size_string

test_bounded_size_string_SOURCES = test_bounded_size_string.cpp


#########################################
# kim::DnaFileReader class test program #
#########################################

check_PROGRAMS += test_dna_file_reader
TESTS          += test_dna_file_reader

test_dna_file_reader_SOURCES = test_dna_file_reader.cpp


########################################
# kim::DnaFileIndex class test program #
########################################

check_PROGRAMS += test_dna_file_index
TESTS          += test_dna_file_index

test_dna_file_index_SOURCES = test_dna_file_index.cpp


######################################
# kim::FileReader class test program #
######################################

check_PROGRAMS += test_file_reader
TESTS          += test_file_reader

test_file_reader_SOURCES = test_file_reader.cpp


#################################
# kim::Hauss class test program #
#################################

check_PROGRAMS += test_gauss
TESTS          += test_gauss

test_gauss_SOURCES = test_gauss.cpp


#####################################
# kim::Exception class test program #
#####################################

check_PROGRAMS += test_kim_exception
TESTS          += test_kim_exception

test_kim_exception_SOURCES = test_kim_exception.cpp


####################################
# kim::Settings class test program #
####################################

check_PROGRAMS += test_kim_settings
TESTS          += test_kim_settings

test_kim_settings_SOURCES = test_kim_settings.cpp


#############################################
# kim::KmerNodesSubindex class test program #
#############################################

check_PROGRAMS += test_kmer_nodes_subindex
TESTS          += test_kmer_nodes_subindex

test_kmer_nodes_subindex_SOURCES = test_kmer_nodes_subindex.cpp


####################################################
# kim::KmerVariantEdgesSubindex class test program #
####################################################

check_PROGRAMS += test_kmer_variant_edges_subindex
TESTS          += test_kmer_variant_edges_subindex

test_kmer_variant_edges_subindex_SOURCES = test_kmer_variant_edges_subindex.cpp


############################################
# kim::KmerVariantGraph class test program #
############################################

check_PROGRAMS += test_kmer_variant_graph
TESTS          += test_kmer_variant_graph

test_kmer_variant_graph_SOURCES = test_kmer_variant_graph.cpp


######################################
# kim::SortHelper class test program #
######################################

check_PROGRAMS += test_sort_helper
TESTS          += test_sort_helper

test_sort_helper_SOURCES = test_sort_helper.cpp

#############################################
# kim::VariantNodesIndex class test program #
#############################################

check_PROGRAMS += test_variant_nodes_index
TESTS          += test_variant_nodes_index

test_variant_nodes_index_SOURCES = test_variant_nodes_index.cpp


##########################################
# kim::VariantFilter* class test program #
##########################################

check_PROGRAMS += test_variant_filter
TESTS          += test_variant_filter

test_variant_filter_SOURCES = test_variant_filter.cpp


#########################
# Dependencies handling #
#########################

# named target
.PHONY: force-rebuild

# rule
force-rebuild:

$(top_builddir)/src/%: force-rebuild
	$(AM_V_at)$(MAKE) $(AM_MAKEFLAGS) "$(@F)" -C "$(@D)"


#################
# Code Coverage #
#################

GCDA_FILES = $(patsubst %.cpp,%.gcda,$(filter %.cpp %.c,$(SOURCES)))
GCNO_FILES = $(GCDA_FILES:.gcda=.gcno)
GCOV_FILES = $(GCDA_FILES:.gcda=.gcov)


####################
# Cleaning targets #
####################

MOSTLYCLEANFILES = *~
CLEANFILES = *~ $(GCDA_FILES) $(GCNO_FILES) $(GCOV_FILES)
DISTCLEANFILES = *~
MAINTAINERCLEANFILES = *~
