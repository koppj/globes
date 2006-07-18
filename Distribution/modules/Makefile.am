INCLUDES = -I$(top_srcdir) -I$(top_srcdir)/lib $(GSLINC) $(LTDLINCL)

if GLB_WO_MODULES
else
pkglib_LTLIBRARIES = glb_prior_module.la
glb_prior_module_la_SOURCES = glb_prior_module.c
glb_prior_module_la_LDFLAGS =  $(top_builddir)/source/libglb.la\
 @LIBS@ -module -avoid-version -no-undefined
endif

lib_LTLIBRARIES = libglobes.la
# Here we try to get a tree like linking scheme including modules
libglobes_la_SOURCES = glb_glb.c glb_modules.c glb_modules.h glb_prior_module.h

if GLB_WO_MODULES
libglobes_la_SOURCES += glb_prior_module.c
endif

libglobes_la_LIBADD = $(top_builddir)/source/libglb.la
libglobes_la_LDFLAGS =  @LIBS@ \
-version-info @GLBLIBVERSION@\
 @LIBLTDL@ 
if GLB_WO_MODULES
else
libglobes_la_LDFLAGS += -dlopen @GLBMODULES@
endif
libglobes_la_CPPFLAGS = -DGLB_MODULE_DIR=\"$(pkglibdir)\"
if GLB_WO_MODULES
else
libglobes_la_DEPENDENCIES = glb_prior_module.la
endif

bin_PROGRAMS = globes
globes_SOURCES =  globes.c 
globes_LDADD = libglobes.la @LIBS@  @LIBLTDL@ -lstdc++ 

if GLB_WO_MODULES
else
globes_LDADD += -dlopen @GLBMODULES@
endif
globes_DEPENDENCIES = libglobes.la
globes_LDADD += $(top_builddir)/lib/libgnu.la

noinst_PROGRAMS = module-example
module_example_SOURCES = module-example.c
module_example_LDADD =  libglobes.la @LIBS@ @LIBLTDL@ -lstdc++
