/* zladiv.f -- translated by f2c (version 20000121).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#ifdef HAVE_G2C_H
#include <g2c.h>
#else
#include <f2c.h>
#endif


/* Double Complex */ VOID zladiv_( ret_val, x, y)
doublecomplex * ret_val;
doublecomplex *x, *y;
{
    /* System generated locals */
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex z__1;

    /* Builtin functions */
    double d_imag();

    /* Local variables */
    static doublereal zi;
    extern /* Subroutine */ int dladiv_();
    static doublereal zr;


/*  -- LAPACK auxiliary routine (version 3.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     October 31, 1992 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  ZLADIV := X / Y, where X and Y are complex.  The computation of X / Y */
/*  will not overflow on an intermediary step unless the results */
/*  overflows. */

/*  Arguments */
/*  ========= */

/*  X       (input) COMPLEX*16 */
/*  Y       (input) COMPLEX*16 */
/*          The complex scalars X and Y. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

    d__1 = x->r;
    d__2 = d_imag(x);
    d__3 = y->r;
    d__4 = d_imag(y);
    dladiv_(&d__1, &d__2, &d__3, &d__4, &zr, &zi);
    z__1.r = zr, z__1.i = zi;
     ret_val->r = z__1.r,  ret_val->i = z__1.i;

    return ;

/*     End of ZLADIV */

} /* zladiv_ */

