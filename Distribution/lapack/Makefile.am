
if WANT_LAPACK
noinst_LTLIBRARIES = libglblapack.la
libglblapack_la_SOURCES = dcabs1.c  dzasum.c  zaxpy.c   zgeev.c   zlacgv.c\
  zlarf.c   zlassq.c  ztrsv.c lsame.c\
dlabad.c  dznrm2.c  zcopy.c   zgehd2.c  zlacpy.c  zlarfb.c  zlatrs.c  zung2r.c\
dladiv.c  idamax.c  zdotc.c   zgehrd.c  zladiv.c  zlarfg.c  zscal.c   zunghr.c\
dlamch.c  ieeeck.c  zdotu.c   zgemm.c   zlahqr.c  zlarft.c  zswap.c   zungqr.c\
dlapy2.c  ilaenv.c  zdscal.c  zgemv.c   zlahrd.c  zlarfx.c  ztrevc.c\
dlapy3.c  izamax.c  zgebak.c  zgerc.c   zlange.c  zlascl.c  ztrmm.c\
dscal.c   xerbla.c  zgebal.c  zhseqr.c  zlanhs.c  zlaset.c  ztrmv.c
libglblapack_la_LDFLAGS = @LIBS@ -static
endif