/* zlarfb.f -- translated by f2c (version 20000121).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#ifdef HAVE_G2C_H
#include <g2c.h>
#else
#include <f2c.h>
#endif


/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;

/* Subroutine */ int zlarfb_(side, trans, direct, storev, m, n, k, v, ldv, t, 
	ldt, c__, ldc, work, ldwork, side_len, trans_len, direct_len, 
	storev_len)
char *side, *trans, *direct, *storev;
integer *m, *n, *k;
doublecomplex *v;
integer *ldv;
doublecomplex *t;
integer *ldt;
doublecomplex *c__;
integer *ldc;
doublecomplex *work;
integer *ldwork;
ftnlen side_len;
ftnlen trans_len;
ftnlen direct_len;
ftnlen storev_len;
{
    /* System generated locals */
    integer c_dim1, c_offset, t_dim1, t_offset, v_dim1, v_offset, work_dim1, 
	    work_offset, i__1, i__2, i__3, i__4, i__5;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    void d_cnjg();

    /* Local variables */
    static integer i__, j;
    extern logical lsame_();
    extern /* Subroutine */ int zgemm_(), zcopy_(), ztrmm_(), zlacgv_();
    static char transt[1];


/*  -- LAPACK auxiliary routine (version 3.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     September 30, 1994 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  ZLARFB applies a complex block reflector H or its transpose H' to a */
/*  complex M-by-N matrix C, from either the left or the right. */

/*  Arguments */
/*  ========= */

/*  SIDE    (input) CHARACTER*1 */
/*          = 'L': apply H or H' from the Left */
/*          = 'R': apply H or H' from the Right */

/*  TRANS   (input) CHARACTER*1 */
/*          = 'N': apply H (No transpose) */
/*          = 'C': apply H' (Conjugate transpose) */

/*  DIRECT  (input) CHARACTER*1 */
/*          Indicates how H is formed from a product of elementary */
/*          reflectors */
/*          = 'F': H = H(1) H(2) . . . H(k) (Forward) */
/*          = 'B': H = H(k) . . . H(2) H(1) (Backward) */

/*  STOREV  (input) CHARACTER*1 */
/*          Indicates how the vectors which define the elementary */
/*          reflectors are stored: */
/*          = 'C': Columnwise */
/*          = 'R': Rowwise */

/*  M       (input) INTEGER */
/*          The number of rows of the matrix C. */

/*  N       (input) INTEGER */
/*          The number of columns of the matrix C. */

/*  K       (input) INTEGER */
/*          The order of the matrix T (= the number of elementary */
/*          reflectors whose product defines the block reflector). */

/*  V       (input) COMPLEX*16 array, dimension */
/*                                (LDV,K) if STOREV = 'C' */
/*                                (LDV,M) if STOREV = 'R' and SIDE = 'L' */
/*                                (LDV,N) if STOREV = 'R' and SIDE = 'R' */
/*          The matrix V. See further details. */

/*  LDV     (input) INTEGER */
/*          The leading dimension of the array V. */
/*          If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M); */
/*          if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N); */
/*          if STOREV = 'R', LDV >= K. */

/*  T       (input) COMPLEX*16 array, dimension (LDT,K) */
/*          The triangular K-by-K matrix T in the representation of the */
/*          block reflector. */

/*  LDT     (input) INTEGER */
/*          The leading dimension of the array T. LDT >= K. */

/*  C       (input/output) COMPLEX*16 array, dimension (LDC,N) */
/*          On entry, the M-by-N matrix C. */
/*          On exit, C is overwritten by H*C or H'*C or C*H or C*H'. */

/*  LDC     (input) INTEGER */
/*          The leading dimension of the array C. LDC >= max(1,M). */

/*  WORK    (workspace) COMPLEX*16 array, dimension (LDWORK,K) */

/*  LDWORK  (input) INTEGER */
/*          The leading dimension of the array WORK. */
/*          If SIDE = 'L', LDWORK >= max(1,N); */
/*          if SIDE = 'R', LDWORK >= max(1,M). */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Quick return if possible */

    /* Parameter adjustments */
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1 * 1;
    v -= v_offset;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1 * 1;
    t -= t_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1 * 1;
    c__ -= c_offset;
    work_dim1 = *ldwork;
    work_offset = 1 + work_dim1 * 1;
    work -= work_offset;

    /* Function Body */
    if (*m <= 0 || *n <= 0) {
	return 0;
    }

    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {
	*(unsigned char *)transt = 'C';
    } else {
	*(unsigned char *)transt = 'N';
    }

    if (lsame_(storev, "C", (ftnlen)1, (ftnlen)1)) {

	if (lsame_(direct, "F", (ftnlen)1, (ftnlen)1)) {

/*           Let  V =  ( V1 )    (first K rows) */
/*                     ( V2 ) */
/*           where  V1  is unit lower triangular. */

	    if (lsame_(side, "L", (ftnlen)1, (ftnlen)1)) {

/*              Form  H * C  or  H' * C  where  C = ( C1 ) */
/*                                                  ( C2 ) */

/*              W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in WORK) */

/*              W := C1' */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    zcopy_(n, &c__[j + c_dim1], ldc, &work[j * work_dim1 + 1],
			     &c__1);
		    zlacgv_(n, &work[j * work_dim1 + 1], &c__1);
/* L10: */
		}

/*              W := W * V1 */

		ztrmm_("Right", "Lower", "No transpose", "Unit", n, k, &c_b1, 
			&v[v_offset], ldv, &work[work_offset], ldwork, (
			ftnlen)5, (ftnlen)5, (ftnlen)12, (ftnlen)4);
		if (*m > *k) {

/*                 W := W + C2'*V2 */

		    i__1 = *m - *k;
		    zgemm_("Conjugate transpose", "No transpose", n, k, &i__1,
			     &c_b1, &c__[*k + 1 + c_dim1], ldc, &v[*k + 1 + 
			    v_dim1], ldv, &c_b1, &work[work_offset], ldwork, (
			    ftnlen)19, (ftnlen)12);
		}

/*              W := W * T'  or  W * T */

		ztrmm_("Right", "Upper", transt, "Non-unit", n, k, &c_b1, &t[
			t_offset], ldt, &work[work_offset], ldwork, (ftnlen)5,
			 (ftnlen)5, (ftnlen)1, (ftnlen)8);

/*              C := C - V * W' */

		if (*m > *k) {

/*                 C2 := C2 - V2 * W' */

		    i__1 = *m - *k;
		    z__1.r = -1., z__1.i = 0.;
		    zgemm_("No transpose", "Conjugate transpose", &i__1, n, k,
			     &z__1, &v[*k + 1 + v_dim1], ldv, &work[
			    work_offset], ldwork, &c_b1, &c__[*k + 1 + c_dim1]
			    , ldc, (ftnlen)12, (ftnlen)19);
		}

/*              W := W * V1' */

		ztrmm_("Right", "Lower", "Conjugate transpose", "Unit", n, k, 
			&c_b1, &v[v_offset], ldv, &work[work_offset], ldwork, 
			(ftnlen)5, (ftnlen)5, (ftnlen)19, (ftnlen)4);

/*              C1 := C1 - W' */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *n;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			i__3 = j + i__ * c_dim1;
			i__4 = j + i__ * c_dim1;
			d_cnjg(&z__2, &work[i__ + j * work_dim1]);
			z__1.r = c__[i__4].r - z__2.r, z__1.i = c__[i__4].i - 
				z__2.i;
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
/* L20: */
		    }
/* L30: */
		}

	    } else if (lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {

/*              Form  C * H  or  C * H'  where  C = ( C1  C2 ) */

/*              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK) */

/*              W := C1 */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    zcopy_(m, &c__[j * c_dim1 + 1], &c__1, &work[j * 
			    work_dim1 + 1], &c__1);
/* L40: */
		}

/*              W := W * V1 */

		ztrmm_("Right", "Lower", "No transpose", "Unit", m, k, &c_b1, 
			&v[v_offset], ldv, &work[work_offset], ldwork, (
			ftnlen)5, (ftnlen)5, (ftnlen)12, (ftnlen)4);
		if (*n > *k) {

/*                 W := W + C2 * V2 */

		    i__1 = *n - *k;
		    zgemm_("No transpose", "No transpose", m, k, &i__1, &c_b1,
			     &c__[(*k + 1) * c_dim1 + 1], ldc, &v[*k + 1 + 
			    v_dim1], ldv, &c_b1, &work[work_offset], ldwork, (
			    ftnlen)12, (ftnlen)12);
		}

/*              W := W * T  or  W * T' */

		ztrmm_("Right", "Upper", trans, "Non-unit", m, k, &c_b1, &t[
			t_offset], ldt, &work[work_offset], ldwork, (ftnlen)5,
			 (ftnlen)5, (ftnlen)1, (ftnlen)8);

/*              C := C - W * V' */

		if (*n > *k) {

/*                 C2 := C2 - W * V2' */

		    i__1 = *n - *k;
		    z__1.r = -1., z__1.i = 0.;
		    zgemm_("No transpose", "Conjugate transpose", m, &i__1, k,
			     &z__1, &work[work_offset], ldwork, &v[*k + 1 + 
			    v_dim1], ldv, &c_b1, &c__[(*k + 1) * c_dim1 + 1], 
			    ldc, (ftnlen)12, (ftnlen)19);
		}

/*              W := W * V1' */

		ztrmm_("Right", "Lower", "Conjugate transpose", "Unit", m, k, 
			&c_b1, &v[v_offset], ldv, &work[work_offset], ldwork, 
			(ftnlen)5, (ftnlen)5, (ftnlen)19, (ftnlen)4);

/*              C1 := C1 - W */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			i__3 = i__ + j * c_dim1;
			i__4 = i__ + j * c_dim1;
			i__5 = i__ + j * work_dim1;
			z__1.r = c__[i__4].r - work[i__5].r, z__1.i = c__[
				i__4].i - work[i__5].i;
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
/* L50: */
		    }
/* L60: */
		}
	    }

	} else {

/*           Let  V =  ( V1 ) */
/*                     ( V2 )    (last K rows) */
/*           where  V2  is unit upper triangular. */

	    if (lsame_(side, "L", (ftnlen)1, (ftnlen)1)) {

/*              Form  H * C  or  H' * C  where  C = ( C1 ) */
/*                                                  ( C2 ) */

/*              W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in WORK) */

/*              W := C2' */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    zcopy_(n, &c__[*m - *k + j + c_dim1], ldc, &work[j * 
			    work_dim1 + 1], &c__1);
		    zlacgv_(n, &work[j * work_dim1 + 1], &c__1);
/* L70: */
		}

/*              W := W * V2 */

		ztrmm_("Right", "Upper", "No transpose", "Unit", n, k, &c_b1, 
			&v[*m - *k + 1 + v_dim1], ldv, &work[work_offset], 
			ldwork, (ftnlen)5, (ftnlen)5, (ftnlen)12, (ftnlen)4);
		if (*m > *k) {

/*                 W := W + C1'*V1 */

		    i__1 = *m - *k;
		    zgemm_("Conjugate transpose", "No transpose", n, k, &i__1,
			     &c_b1, &c__[c_offset], ldc, &v[v_offset], ldv, &
			    c_b1, &work[work_offset], ldwork, (ftnlen)19, (
			    ftnlen)12);
		}

/*              W := W * T'  or  W * T */

		ztrmm_("Right", "Lower", transt, "Non-unit", n, k, &c_b1, &t[
			t_offset], ldt, &work[work_offset], ldwork, (ftnlen)5,
			 (ftnlen)5, (ftnlen)1, (ftnlen)8);

/*              C := C - V * W' */

		if (*m > *k) {

/*                 C1 := C1 - V1 * W' */

		    i__1 = *m - *k;
		    z__1.r = -1., z__1.i = 0.;
		    zgemm_("No transpose", "Conjugate transpose", &i__1, n, k,
			     &z__1, &v[v_offset], ldv, &work[work_offset], 
			    ldwork, &c_b1, &c__[c_offset], ldc, (ftnlen)12, (
			    ftnlen)19);
		}

/*              W := W * V2' */

		ztrmm_("Right", "Upper", "Conjugate transpose", "Unit", n, k, 
			&c_b1, &v[*m - *k + 1 + v_dim1], ldv, &work[
			work_offset], ldwork, (ftnlen)5, (ftnlen)5, (ftnlen)
			19, (ftnlen)4);

/*              C2 := C2 - W' */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *n;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			i__3 = *m - *k + j + i__ * c_dim1;
			i__4 = *m - *k + j + i__ * c_dim1;
			d_cnjg(&z__2, &work[i__ + j * work_dim1]);
			z__1.r = c__[i__4].r - z__2.r, z__1.i = c__[i__4].i - 
				z__2.i;
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
/* L80: */
		    }
/* L90: */
		}

	    } else if (lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {

/*              Form  C * H  or  C * H'  where  C = ( C1  C2 ) */

/*              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK) */

/*              W := C2 */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    zcopy_(m, &c__[(*n - *k + j) * c_dim1 + 1], &c__1, &work[
			    j * work_dim1 + 1], &c__1);
/* L100: */
		}

/*              W := W * V2 */

		ztrmm_("Right", "Upper", "No transpose", "Unit", m, k, &c_b1, 
			&v[*n - *k + 1 + v_dim1], ldv, &work[work_offset], 
			ldwork, (ftnlen)5, (ftnlen)5, (ftnlen)12, (ftnlen)4);
		if (*n > *k) {

/*                 W := W + C1 * V1 */

		    i__1 = *n - *k;
		    zgemm_("No transpose", "No transpose", m, k, &i__1, &c_b1,
			     &c__[c_offset], ldc, &v[v_offset], ldv, &c_b1, &
			    work[work_offset], ldwork, (ftnlen)12, (ftnlen)12)
			    ;
		}

/*              W := W * T  or  W * T' */

		ztrmm_("Right", "Lower", trans, "Non-unit", m, k, &c_b1, &t[
			t_offset], ldt, &work[work_offset], ldwork, (ftnlen)5,
			 (ftnlen)5, (ftnlen)1, (ftnlen)8);

/*              C := C - W * V' */

		if (*n > *k) {

/*                 C1 := C1 - W * V1' */

		    i__1 = *n - *k;
		    z__1.r = -1., z__1.i = 0.;
		    zgemm_("No transpose", "Conjugate transpose", m, &i__1, k,
			     &z__1, &work[work_offset], ldwork, &v[v_offset], 
			    ldv, &c_b1, &c__[c_offset], ldc, (ftnlen)12, (
			    ftnlen)19);
		}

/*              W := W * V2' */

		ztrmm_("Right", "Upper", "Conjugate transpose", "Unit", m, k, 
			&c_b1, &v[*n - *k + 1 + v_dim1], ldv, &work[
			work_offset], ldwork, (ftnlen)5, (ftnlen)5, (ftnlen)
			19, (ftnlen)4);

/*              C2 := C2 - W */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			i__3 = i__ + (*n - *k + j) * c_dim1;
			i__4 = i__ + (*n - *k + j) * c_dim1;
			i__5 = i__ + j * work_dim1;
			z__1.r = c__[i__4].r - work[i__5].r, z__1.i = c__[
				i__4].i - work[i__5].i;
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
/* L110: */
		    }
/* L120: */
		}
	    }
	}

    } else if (lsame_(storev, "R", (ftnlen)1, (ftnlen)1)) {

	if (lsame_(direct, "F", (ftnlen)1, (ftnlen)1)) {

/*           Let  V =  ( V1  V2 )    (V1: first K columns) */
/*           where  V1  is unit upper triangular. */

	    if (lsame_(side, "L", (ftnlen)1, (ftnlen)1)) {

/*              Form  H * C  or  H' * C  where  C = ( C1 ) */
/*                                                  ( C2 ) */

/*              W := C' * V'  =  (C1'*V1' + C2'*V2') (stored in WORK) */

/*              W := C1' */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    zcopy_(n, &c__[j + c_dim1], ldc, &work[j * work_dim1 + 1],
			     &c__1);
		    zlacgv_(n, &work[j * work_dim1 + 1], &c__1);
/* L130: */
		}

/*              W := W * V1' */

		ztrmm_("Right", "Upper", "Conjugate transpose", "Unit", n, k, 
			&c_b1, &v[v_offset], ldv, &work[work_offset], ldwork, 
			(ftnlen)5, (ftnlen)5, (ftnlen)19, (ftnlen)4);
		if (*m > *k) {

/*                 W := W + C2'*V2' */

		    i__1 = *m - *k;
		    zgemm_("Conjugate transpose", "Conjugate transpose", n, k,
			     &i__1, &c_b1, &c__[*k + 1 + c_dim1], ldc, &v[(*k 
			    + 1) * v_dim1 + 1], ldv, &c_b1, &work[work_offset]
			    , ldwork, (ftnlen)19, (ftnlen)19);
		}

/*              W := W * T'  or  W * T */

		ztrmm_("Right", "Upper", transt, "Non-unit", n, k, &c_b1, &t[
			t_offset], ldt, &work[work_offset], ldwork, (ftnlen)5,
			 (ftnlen)5, (ftnlen)1, (ftnlen)8);

/*              C := C - V' * W' */

		if (*m > *k) {

/*                 C2 := C2 - V2' * W' */

		    i__1 = *m - *k;
		    z__1.r = -1., z__1.i = 0.;
		    zgemm_("Conjugate transpose", "Conjugate transpose", &
			    i__1, n, k, &z__1, &v[(*k + 1) * v_dim1 + 1], ldv,
			     &work[work_offset], ldwork, &c_b1, &c__[*k + 1 + 
			    c_dim1], ldc, (ftnlen)19, (ftnlen)19);
		}

/*              W := W * V1 */

		ztrmm_("Right", "Upper", "No transpose", "Unit", n, k, &c_b1, 
			&v[v_offset], ldv, &work[work_offset], ldwork, (
			ftnlen)5, (ftnlen)5, (ftnlen)12, (ftnlen)4);

/*              C1 := C1 - W' */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *n;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			i__3 = j + i__ * c_dim1;
			i__4 = j + i__ * c_dim1;
			d_cnjg(&z__2, &work[i__ + j * work_dim1]);
			z__1.r = c__[i__4].r - z__2.r, z__1.i = c__[i__4].i - 
				z__2.i;
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
/* L140: */
		    }
/* L150: */
		}

	    } else if (lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {

/*              Form  C * H  or  C * H'  where  C = ( C1  C2 ) */

/*              W := C * V'  =  (C1*V1' + C2*V2')  (stored in WORK) */

/*              W := C1 */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    zcopy_(m, &c__[j * c_dim1 + 1], &c__1, &work[j * 
			    work_dim1 + 1], &c__1);
/* L160: */
		}

/*              W := W * V1' */

		ztrmm_("Right", "Upper", "Conjugate transpose", "Unit", m, k, 
			&c_b1, &v[v_offset], ldv, &work[work_offset], ldwork, 
			(ftnlen)5, (ftnlen)5, (ftnlen)19, (ftnlen)4);
		if (*n > *k) {

/*                 W := W + C2 * V2' */

		    i__1 = *n - *k;
		    zgemm_("No transpose", "Conjugate transpose", m, k, &i__1,
			     &c_b1, &c__[(*k + 1) * c_dim1 + 1], ldc, &v[(*k 
			    + 1) * v_dim1 + 1], ldv, &c_b1, &work[work_offset]
			    , ldwork, (ftnlen)12, (ftnlen)19);
		}

/*              W := W * T  or  W * T' */

		ztrmm_("Right", "Upper", trans, "Non-unit", m, k, &c_b1, &t[
			t_offset], ldt, &work[work_offset], ldwork, (ftnlen)5,
			 (ftnlen)5, (ftnlen)1, (ftnlen)8);

/*              C := C - W * V */

		if (*n > *k) {

/*                 C2 := C2 - W * V2 */

		    i__1 = *n - *k;
		    z__1.r = -1., z__1.i = 0.;
		    zgemm_("No transpose", "No transpose", m, &i__1, k, &z__1,
			     &work[work_offset], ldwork, &v[(*k + 1) * v_dim1 
			    + 1], ldv, &c_b1, &c__[(*k + 1) * c_dim1 + 1], 
			    ldc, (ftnlen)12, (ftnlen)12);
		}

/*              W := W * V1 */

		ztrmm_("Right", "Upper", "No transpose", "Unit", m, k, &c_b1, 
			&v[v_offset], ldv, &work[work_offset], ldwork, (
			ftnlen)5, (ftnlen)5, (ftnlen)12, (ftnlen)4);

/*              C1 := C1 - W */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			i__3 = i__ + j * c_dim1;
			i__4 = i__ + j * c_dim1;
			i__5 = i__ + j * work_dim1;
			z__1.r = c__[i__4].r - work[i__5].r, z__1.i = c__[
				i__4].i - work[i__5].i;
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
/* L170: */
		    }
/* L180: */
		}

	    }

	} else {

/*           Let  V =  ( V1  V2 )    (V2: last K columns) */
/*           where  V2  is unit lower triangular. */

	    if (lsame_(side, "L", (ftnlen)1, (ftnlen)1)) {

/*              Form  H * C  or  H' * C  where  C = ( C1 ) */
/*                                                  ( C2 ) */

/*              W := C' * V'  =  (C1'*V1' + C2'*V2') (stored in WORK) */

/*              W := C2' */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    zcopy_(n, &c__[*m - *k + j + c_dim1], ldc, &work[j * 
			    work_dim1 + 1], &c__1);
		    zlacgv_(n, &work[j * work_dim1 + 1], &c__1);
/* L190: */
		}

/*              W := W * V2' */

		ztrmm_("Right", "Lower", "Conjugate transpose", "Unit", n, k, 
			&c_b1, &v[(*m - *k + 1) * v_dim1 + 1], ldv, &work[
			work_offset], ldwork, (ftnlen)5, (ftnlen)5, (ftnlen)
			19, (ftnlen)4);
		if (*m > *k) {

/*                 W := W + C1'*V1' */

		    i__1 = *m - *k;
		    zgemm_("Conjugate transpose", "Conjugate transpose", n, k,
			     &i__1, &c_b1, &c__[c_offset], ldc, &v[v_offset], 
			    ldv, &c_b1, &work[work_offset], ldwork, (ftnlen)
			    19, (ftnlen)19);
		}

/*              W := W * T'  or  W * T */

		ztrmm_("Right", "Lower", transt, "Non-unit", n, k, &c_b1, &t[
			t_offset], ldt, &work[work_offset], ldwork, (ftnlen)5,
			 (ftnlen)5, (ftnlen)1, (ftnlen)8);

/*              C := C - V' * W' */

		if (*m > *k) {

/*                 C1 := C1 - V1' * W' */

		    i__1 = *m - *k;
		    z__1.r = -1., z__1.i = 0.;
		    zgemm_("Conjugate transpose", "Conjugate transpose", &
			    i__1, n, k, &z__1, &v[v_offset], ldv, &work[
			    work_offset], ldwork, &c_b1, &c__[c_offset], ldc, 
			    (ftnlen)19, (ftnlen)19);
		}

/*              W := W * V2 */

		ztrmm_("Right", "Lower", "No transpose", "Unit", n, k, &c_b1, 
			&v[(*m - *k + 1) * v_dim1 + 1], ldv, &work[
			work_offset], ldwork, (ftnlen)5, (ftnlen)5, (ftnlen)
			12, (ftnlen)4);

/*              C2 := C2 - W' */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *n;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			i__3 = *m - *k + j + i__ * c_dim1;
			i__4 = *m - *k + j + i__ * c_dim1;
			d_cnjg(&z__2, &work[i__ + j * work_dim1]);
			z__1.r = c__[i__4].r - z__2.r, z__1.i = c__[i__4].i - 
				z__2.i;
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
/* L200: */
		    }
/* L210: */
		}

	    } else if (lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {

/*              Form  C * H  or  C * H'  where  C = ( C1  C2 ) */

/*              W := C * V'  =  (C1*V1' + C2*V2')  (stored in WORK) */

/*              W := C2 */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    zcopy_(m, &c__[(*n - *k + j) * c_dim1 + 1], &c__1, &work[
			    j * work_dim1 + 1], &c__1);
/* L220: */
		}

/*              W := W * V2' */

		ztrmm_("Right", "Lower", "Conjugate transpose", "Unit", m, k, 
			&c_b1, &v[(*n - *k + 1) * v_dim1 + 1], ldv, &work[
			work_offset], ldwork, (ftnlen)5, (ftnlen)5, (ftnlen)
			19, (ftnlen)4);
		if (*n > *k) {

/*                 W := W + C1 * V1' */

		    i__1 = *n - *k;
		    zgemm_("No transpose", "Conjugate transpose", m, k, &i__1,
			     &c_b1, &c__[c_offset], ldc, &v[v_offset], ldv, &
			    c_b1, &work[work_offset], ldwork, (ftnlen)12, (
			    ftnlen)19);
		}

/*              W := W * T  or  W * T' */

		ztrmm_("Right", "Lower", trans, "Non-unit", m, k, &c_b1, &t[
			t_offset], ldt, &work[work_offset], ldwork, (ftnlen)5,
			 (ftnlen)5, (ftnlen)1, (ftnlen)8);

/*              C := C - W * V */

		if (*n > *k) {

/*                 C1 := C1 - W * V1 */

		    i__1 = *n - *k;
		    z__1.r = -1., z__1.i = 0.;
		    zgemm_("No transpose", "No transpose", m, &i__1, k, &z__1,
			     &work[work_offset], ldwork, &v[v_offset], ldv, &
			    c_b1, &c__[c_offset], ldc, (ftnlen)12, (ftnlen)12)
			    ;
		}

/*              W := W * V2 */

		ztrmm_("Right", "Lower", "No transpose", "Unit", m, k, &c_b1, 
			&v[(*n - *k + 1) * v_dim1 + 1], ldv, &work[
			work_offset], ldwork, (ftnlen)5, (ftnlen)5, (ftnlen)
			12, (ftnlen)4);

/*              C1 := C1 - W */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			i__3 = i__ + (*n - *k + j) * c_dim1;
			i__4 = i__ + (*n - *k + j) * c_dim1;
			i__5 = i__ + j * work_dim1;
			z__1.r = c__[i__4].r - work[i__5].r, z__1.i = c__[
				i__4].i - work[i__5].i;
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
/* L230: */
		    }
/* L240: */
		}

	    }

	}
    }

    return 0;

/*     End of ZLARFB */

} /* zlarfb_ */
