## Process this file with automake to produce Makefile.in

INCLUDES=-I$(top_builddir)

AUTOMAKE_OPTIONS = no-dependencies foreign

if INSTALL_LTDL
include_HEADERS = ltdl.h
lib_LTLIBRARIES = libltdl.la
else
noinst_HEADERS = ltdl.h
endif

if CONVENIENCE_LTDL
noinst_LTLIBRARIES = libltdlc.la
endif

## Make sure these will be cleaned even when they're not built by
## default.
CLEANFILES = libltdl.la libltdlc.la

libltdl_la_SOURCES = ltdl.c
libltdl_la_LDFLAGS = -no-undefined -version-info 4:0:1
libltdl_la_LIBADD = $(LIBADD_DL)

libltdlc_la_SOURCES = ltdl.c
libltdlc_la_LIBADD = $(LIBADD_DL)
