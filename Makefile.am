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

SUBDIRS = resources src tests

if BUILD_DOC

  SUBDIRS += doc

if MAKE_RELEASE

  SUBDIRS += packaging

# Fake target to force packages construction
always-build:

else

# Fake target to abort packages construction
always-build:
	@echo
	@echo "====================================================================="
	@echo
	@echo "Release construction is disabled. Your local system seems to miss"
	@echo "the 'dpkg' command."
	@echo
	@echo "====================================================================="
	@echo
	@false

endif

else

# Fake target to abort packages construction
always-build:
	@echo
	@echo "====================================================================="
	@echo
	@echo "Release construction is disabled. Your local system probably misses"
	@echo "one of the 'doxygen' or 'pdftex' command."
	@echo
	@echo "====================================================================="
	@echo
	@false

endif

dist_pkgdata_DATA = \
  AUTHORS ChangeLog README.md \
  LICENSE.md LICENSE-fr.md LICENSE-en.md

nobase_include_HEADERS = @top_srcdir@/BitMagic/src/*.h

EXTRA_DIST = \
  COPYING INSTALL NEWS \
  config/git-version-gen config/version config/version.configure.ac

ACLOCAL_AMFLAGS = -I config/m4

deb: @PACKAGE@_@PACKAGE_VERSION@_@ARCH@.deb

@PACKAGE@_@PACKAGE_VERSION@_@ARCH@.deb: always-build
	@if [ "@DEB_BUILDER@" == "" ]; then \
	  echo "You need to install a debian package builder and to reconfigure this package.";\
	  exit 1;\
	fi
	make dist
	$(am__remove_distdir)
	GZIP=$(GZIP_ENV) gunzip -c $(distdir).tar.gz | $(AMTAR) xf -
	mkdir $(distdir)/_build
	mkdir -p $(distdir)/_inst/_deb/$(prefix)
	chmod a-w $(distdir)
	dc_install_base=`$(am__cd) $(distdir)/_inst/_deb/$(prefix) && pwd | sed -e 's,^[^:\\/]:[\\/],/,'` \
	  && dc_destdir="$${TMPDIR-/tmp}/am-dc-$$$$/" \
	  && cd $(distdir)/_build \
	  && HACKED_PREFIX=$(prefix) ../configure --srcdir=.. --prefix="$$dc_install_base" CC="$(CC)" CXX="$(CXX)" \
	    $(DISTCHECK_CONFIGURE_FLAGS) \
	  && $(MAKE) $(AM_MAKEFLAGS) localedir=$$localedir \
	  && $(MAKE) $(AM_MAKEFLAGS) check \
	  && $(MAKE) $(AM_MAKEFLAGS) install \
	  && $(MAKE) $(AM_MAKEFLAGS) installcheck \
	  && cd ../_inst/_deb \
	  && cp "$(abs_top_srcdir)/resources/@PACKAGE@-@PACKAGE_VERSION@.stats" "./$(pkgdatadir)/" \
	  && find . -type f -exec md5sum {} \; | sed "s,\.\/,,g" > ../md5sums \
          && mkdir DEBIAN \
	  && mv ../md5sums DEBIAN/ \
	  && cp ../../_build/packaging/postinst DEBIAN/ \
	  && chmod 755 DEBIAN/postinst \
	  && cp ../../_build/packaging/prerm DEBIAN/ \
	  && chmod 755 DEBIAN/prerm \
	  && sed "s/@INSTALLED_SIZE@/`du -ck . | tail -1 | cut -f1`/" \
	       ../../_build/packaging/control > DEBIAN/control \
	  && cd ../ \
	  && dpkg -b _deb $@ \
	  && mv $@ ../../ \
	  && rm -f "$$dc_install_base/DEBIAN/control" \
	  && rm -f "$$dc_install_base/DEBIAN/prerm" \
	  && rm -f "$$dc_install_base/DEBIAN/postinst" \
	  && cd ../_build \
	  && $(MAKE) $(AM_MAKEFLAGS) uninstall \
	  && $(MAKE) $(AM_MAKEFLAGS) distuninstallcheck_dir="$$dc_install_base" \
	        distuninstallcheck \
	  && chmod -R a-w "$$dc_install_base" \
	  && ({ \
	       (cd ../.. && $(mkinstalldirs) "$$dc_destdir") \
	       && $(MAKE) $(AM_MAKEFLAGS) DESTDIR="$$dc_destdir" install \
	       && $(MAKE) $(AM_MAKEFLAGS) DESTDIR="$$dc_destdir" uninstall \
	       && $(MAKE) $(AM_MAKEFLAGS) DESTDIR="$$dc_destdir" \
	            distuninstallcheck_dir="$$dc_destdir" distuninstallcheck; \
	      } \
	      || { rm -rf "$$dc_destdir"; exit 1; }) \
	  && rm -rf "$$dc_destdir" \
	  && $(MAKE) $(AM_MAKEFLAGS) distcleancheck
	$(am__remove_distdir)
	@echo "$@ is ready for distribution" | \
	  sed 'h;s/./=/g;p;x;p;x'


MOSTLYCLEANFILES=*~ config/*~
CLEANFILES=*~ config/*~
DISTCLEANFILES=*~ $(distdir).tar.gz
DISTCLEANFILES+=@PACKAGE@_@PACKAGE_VERSION@_@ARCH@.deb
MAINTAINERCLEANFILES=*~ config/*~
