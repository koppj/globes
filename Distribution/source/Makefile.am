INCLUDES = -I$(top_srcdir)

AM_YFLAGS =  -d 


lib_LTLIBRARIES = libglobes.la
libglobes_la_SOURCES =  glb_minimize.c glb_error.c glb_path.c glb_multiex.c\
 glb_version.c \
glb_rate_engine.c  \
  glb_parser.y  glb_min_sup.c glb_fluxes.c glb_smear.c glb_lexer.l glb_probability.cc glb_profile.c\
  user_chi.c glb_wrapper.c glb_modules.c\
glb_minimize.h  glb_error.h glb_path.h glb_parser_addons.h glb_multiex.h glb_probability.h glb_smear.h \
glb_rate_engine.h glb_types.h glb_parser.h glb_min_sup.h  glb_version.h \
glb_lexer.h glb_fluxes.h glb_parser_type.h glb_wrapper.h glb_modules.h\
glb_profile.h


# libglobes_la_LIBADD =  
libglobes_la_LDFLAGS =  @LIBS@ -version-info @GLBLIBVERSION@
libglobes_la_CFLAGS= $(AM_CFLAGS)
libglobes_la_CXXFLAGS= $(AM_CFLAGS)


pkginclude_HEADERS = ../globes/globes.h

bin_PROGRAMS = globes
globes_SOURCES = glb_min_sup.c glb_minimize.c glb_error.c glb_path.c\
glb_multiex.c \
glb_version.c \
glb_rate_engine.c glb_parser.y glb_fluxes.c glb_smear.c glb_lexer.l glb_probability.cc glb_profile.c\
 user_chi.c \
glb_wrapper.c glb_modules.c globes.c 
globes_CXXFLAGS= $(AM_CFLAGS)
globes_CFLAGS= $(AM_CFLAGS)
globes_LDADD = @LIBS@ 
