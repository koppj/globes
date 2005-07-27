ACLOCAL_AMFLAGS = -I m4

DIST_SUBDIRS = lapack source data libltdl modules globes libf2c
SUBDIRS = libf2c lapack libltdl source data modules globes 
INCLUDES = $(LTDLINCL) $(F2CINC)

EXTRA_DIST = examples/Makefile.in examples/example1.c\
  examples/example3.c examples/myio.c\
examples/example-tour.c  examples/example2.c  examples/example4.c\
examples/myio.h examples/XCC.dat examples/NuFact.glb\
BUGS globes-script lapack/Notice libf2c/Notice

nodist_bin_SCRIPTS = globes-config globes-module
noinst_SCRIPTS = globes-script


globes-config: globes-config.src
	$(top_srcdir)/globes-script globes-config.src globes-config

globes-module: globes-module.src
	$(top_srcdir)/globes-script globes-module.src globes-module 

CLEANFILES = globes-module globes-config

AM_YFLAGS =  -d 

if MAKE_RPMS

rpm: $(DISTFILES) 
	${MAKE} dist
# adding my own helpless stuff ...
	cp $(RPM_TARBALL) $(RPM_SOURCE_DIR)/
	$(RPM_PROG) $(RPM_ARGS) --clean @PACKAGE@.spec
	$(RPM_PROG) --rmsource @PACKAGE@.spec
	cp $(RPM_SRPM_NAME) ./
	cp $(RPM_TARGET) ./
	rm -f $(RPM_SRPM_NAME)
	rm -f $(RPM_SOURCE_DIR)/$(RPM_TARBALL)
	rm -rf $(RPM_BUILD_DIR)/$(RPM_VERSION)
	rm -f $(RPM_TARGET)
	@(echo "$(distdir) rpms ready for distribution: "; \
          list='$(RPM_CLEAN_NAMES)'; for i in $$list; do echo $$i; done) | \
          sed -e '1{h;s/./=/g;p;x;}' -e '$${p;x;}'

rpm-clean:
	rm -f $(RPM_SRPM_NAME)
	rm -f $(RPM_SOURCE_DIR)/$(RPM_TARBALL)
	rm -rf $(RPM_BUILD_DIR)/$(RPM_VERSION)
	rm -f $(RPM_TARGET)

.PHONY: rpm rpm-clean

else
endif
