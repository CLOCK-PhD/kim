@DX_RULES@

EXTRA_DATA = \
	$(DX_CLEAN_HTML) \
	$(DX_CLEAN_CHM) \
	$(DX_CLEAN_CHI) \
	$(DX_CLEAN_RTF) \
	$(DX_CLEAN_XML) \
	$(DX_CLEAN_PS) \
	$(DX_CLEAN_PDF)

if BUILD_DOC_MAN

man3_MANS = @DX_DOCDIR@/man/man3/gkampi_*

$(man3_MANS): doxygen-doc

endif

#### HTML ####

install-data-hook: $(EXTRA_DATA)
	@$(NORMAL_INSTALL)
	$(AM_V_at)test -n "$(docdir)" || exit 1
	$(AM_V_at)for p in $(EXTRA_DATA); do \
	  if test -f "$$p"; then \
	    dirname "$$p"; \
	  else \
	    find "$$p" -type d -print; \
	  fi; \
	done | $(SED) 's,^$(DX_DOCDIR),,' | \
	while read dest; do \
	  echo " $(MKDIR_P) '$(DESTDIR)$(docdir)/$$dest'"; \
	  $(MKDIR_P) "$(DESTDIR)$(docdir)/$$dest" || exit 1; \
	done
	$(AM_V_at)list='$(EXTRA_DATA)'; \
	for p in $$list; do \
	  if test -f "$$p"; then \
	    echo "$$p"; \
	  else \
	    find "$$p" -type f -print; \
	  fi; \
	done | \
	while read file; do \
	  dest=$(DESTDIR)$(docdir)/$$(dirname "$${file#$(DX_DOCDIR)/}"); \
	  echo " $(INSTALL_DATA) '$$file' '$$dest'"; \
	  $(INSTALL_DATA) "$$file" "$$dest" || exit $$?; \
	done

uninstall-hook:
	@$(NORMAL_UNINSTALL)
	$(AM_V_at)list='$(EXTRA_DATA)'; test -n "$(docdir)" || list=; \
	files=`for p in $$list; do echo $$p; done | sed -e 's|^.*/||'`; \
	dir='$(DESTDIR)$(docdir)'; \
	test -z "$$files" \
	|| { test ! -d "$$dir" && test ! -f "$$dir" && test ! -r "$$dir"; } \
	|| { echo " ( cd '$$dir' && rm -rf" $$files ")"; \
	         $(am__cd) "$$dir" && rm -rf $$files; }

all-local: doxygen-doc

MOSTLYCLEANFILES=*~ $(DX_CLEANFILES)
CLEANFILES=*~
DISTCLEANFILES=*~
MAINTAINERCLEANFILES=*~
