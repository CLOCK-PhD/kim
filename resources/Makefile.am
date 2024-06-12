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

STATS_FILE_BASENAME = @PACKAGE@-@PACKAGE_VERSION@
CLOC_FILE = $(STATS_FILE_BASENAME).cloc
OHCOUNT_FILE = $(STATS_FILE_BASENAME).ohcount
LIST_FILE = $(STATS_FILE_BASENAME).filelist
STATS_FILE = $(STATS_FILE_BASENAME).stats

dist_pkgdata_DATA =	\
  kim.bash_completion	\
  test-reads.fastq	\
  test-sequences.fasta	\
  $(STATS_FILE)

$(LIST_FILE):
	@echo "Generate package file list"
	$(AM_V_at)git ls-tree --name-only -r  HEAD "$(abs_top_srcdir)" > "$@" || true

if CLOC

$(CLOC_FILE): $(LIST_FILE)
	@echo "Generate cloc summary"
	$(AM_V_at)echo "CLoC summary" | sed 'h;s/./=/g;p;x;p;x' > "$@"
	$(AM_V_at)echo >> "$@"
	$(AM_V_at)@CLOC@ --list-file="$<" --quiet | tail -n +2 >> "$@"

else

$(CLOC_FILE):
	@echo "You should install cloc (source code line counter) and reconfigure this package."
	$(AM_V_at)echo "CLoC isn't installed. No report created" > "$@"

endif

if OHCOUNT

$(OHCOUNT_FILE): $(LIST_FILE)
	@echo "Generate OhCount summary"
	$(AM_V_at)echo "OhCount summary" | sed 'h;s/./=/g;p;x;p;x' > "$@"
	$(AM_V_at)@OHCOUNT@ $$(cat "$<") | tail -n +2 >> "$@"

else

$(OHCOUNT_FILE): $(LIST_FILE)
	@echo "You should install ohcount (source code line counter) and reconfigure this package."
	$(AM_V_at)echo "OhCount isn't installed. No report created" > "$@"

endif

force-renew:

$(STATS_FILE): $(CLOC_FILE) $(OHCOUNT_FILE) $(LIST_FILE) force-renew
	@echo "Backup existing code statistics (if any)"
	$(AM_V_at)test ! -f "$(abs_srcdir)/$(STATS_FILE)" || cat "$(abs_srcdir)/$(STATS_FILE)" > "$(STATS_FILE).tmp"
	@echo "Generate code statistics"
	$(AM_V_at)echo "# $@ #" | sed 'h;s/./#/g;p;x;p;x' > "$@"
	$(AM_V_at)echo >> "$@"
	$(AM_V_at)cat $(CLOC_FILE) >> "$@"
	$(AM_V_at)echo >> "$@"
	$(AM_V_at)cat $(OHCOUNT_FILE) >> "$@"
	$(AM_V_at)echo >> "$@"
	$(AM_V_at)printf "%.78d\n" | tr 0 - >> "$@"
	$(AM_V_at)echo >> "$@"
	$(AM_V_at)if test ! -s "$(LIST_FILE)" -a -f "$(STATS_FILE).tmp"; then		\
	  echo "Using existing code statistics backup file since generation failed";	\
	  mv "$(STATS_FILE).tmp" "$(STATS_FILE)";					\
	fi
	@echo "Removing temporary files and backup (if any)"
	$(AM_V_at)$(RM) "$(STATS_FILE).tmp" "$(CLOC_FILE)" "$(OHCOUNT_FILE)" "$(LIST_FILE)"
	$(AM_V_at)cat "$(STATS_FILE)"

cloc: $(CLOC_FILE)

ohcount: $(OHCOUNT_FILE)

filelist: $(LIST_FILE)

stats: $(STATS_FILE)

MOSTLYCLEANFILES=*~
CLEANFILES=*~ $(LIST_FILE) $(CLOC_FILE) $(OHCOUNT_FILE) $(STATS_FILE)
DISTCLEANFILES=*~
MAINTAINERCLEANFILES=*~
