/* dcabs1.f -- translated by f2c (version 20000121).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#ifdef HAVE_G2C_H
#include <g2c.h>
#else
#include <f2c.h>
#endif


doublereal dcabs1_(z__)
doublecomplex *z__;
{
    /* System generated locals */
    doublereal ret_val;
    static doublecomplex equiv_0[1];

    /* Local variables */
#define t ((doublereal *)equiv_0)
#define zz (equiv_0)

    zz->r = z__->r, zz->i = z__->i;
    ret_val = abs(t[0]) + abs(t[1]);
    return ret_val;
} /* dcabs1_ */

#undef zz
#undef t

