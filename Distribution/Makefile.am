DIST_SUBDIRS = lapack source data libltdl modules globes libf2c
SUBDIRS = libf2c lapack libltdl source data modules globes 
INCLUDES = $(LTDLINCL) $(F2CINCL)

EXTRA_DIST = examples/Makefile.in examples/example1.c\
  examples/example3.c examples/myio.c\
examples/example-tour.c  examples/example2.c  examples/example4.c\
examples/myio.h examples/XCC.dat examples/NuFact.glb\
BUGS globes-script lapack/COPYING

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
	@echo =================================
	@echo Congratulations, $(RPM_TARGET) "(and friends)" should now exist.
	@echo =================================


rpm-clean:
	rm -f $(RPM_SOURCE_DIR)/$(RPM_TARBALL)
	rm -rf $(RPM_BUILD_DIR)/$(RPM_VERSION)
	rm -f $(RPM_TARGET)

.PHONY: rpm rpm-clean

else
endif
