SUBDIRS = src

if BUILD_DOC

  SUBDIRS += doc

endif

dist_pkgdata_DATA = \
  AUTHORS ChangeLog README.md \
  LICENSE-CeCILL-en.md LICENSE-CeCILL-fr.md \
  LICENSE.md

EXTRA_DIST = \
  COPYING INSTALL NEWS

ACLOCAL_AMFLAGS = -I config/m4

MOSTLYCLEANFILES=*~ config/*~
CLEANFILES=*~ config/*~
DISTCLEANFILES=*~ $(distdir).tar.gz
MAINTAINERCLEANFILES=*~ config/*~
