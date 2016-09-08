/* erf64.f -- translated by f2c (version 20090411).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"


/* $Id: erf64.F,v 1.1.1.1 1996/04/01 15:01:53 mclareni Exp $ */

/* $Log: erf64.F,v $ */
/* Revision 1.1.1.1  1996/04/01 15:01:53  mclareni */
/* Mathlib gen */


doublereal derf_0_(int n__, doublereal *x)
{
    /* Initialized data */

    static doublereal p1[4] = { 242.66795523053175,21.979261618294152,
	    6.9963834886191355,-.035609843701815385 };
    static doublereal q1[4] = { 215.0588758698612,91.164905404514901,
	    15.082797630407787,1. };
    static doublereal p2[8] = { 300.459261020161601,451.918953711872942,
	    339.320816734343687,152.989285046940404,43.1622272220567353,
	    7.21175825088309366,.564195517478973971,-1.36864857382716707e-7 };
    static doublereal q2[8] = { 300.459260956983293,790.950925327898027,
	    931.354094850609621,638.980264465631167,277.585444743987643,
	    77.000152935229473,12.7827273196294235,1. };
    static doublereal p3[5] = { -.00299610707703542174,-.0494730910623250734,
	    -.22695659353968693,-.278661308609647788,-.0223192459734184686 };
    static doublereal q3[5] = { .0106209230528467918,.191308926107829841,
	    1.05167510706793207,1.98733201817135256,1. };

    /* System generated locals */
    doublereal ret_val, d__1;

    /* Builtin functions */
    double exp(doublereal);

    /* Local variables */
    static doublereal h__;
    static integer i__;
    static doublereal v, y, hc, ap, aq;
    static logical lef;



/* $Id: imp64.inc,v 1.1.1.1 1996/04/01 15:02:59 mclareni Exp $ */

/* $Log: imp64.inc,v $ */
/* Revision 1.1.1.1  1996/04/01 15:02:59  mclareni */
/* Mathlib gen */


/* imp64.inc */



    switch(n__) {
	case 1: goto L_derfc;
	}

    lef = TRUE_;
    goto L9;

L_derfc:
    lef = FALSE_;
L9:
    v = abs(*x);
    if (v < .5) {
/* Computing 2nd power */
	d__1 = v;
	y = d__1 * d__1;
	ap = p1[3];
	aq = q1[3];
	for (i__ = 2; i__ >= 0; --i__) {
	    ap = p1[i__] + y * ap;
/* L1: */
	    aq = q1[i__] + y * aq;
	}
	h__ = *x * ap / aq;
	hc = 1 - h__;
    } else {
	if (v < 4.) {
	    ap = p2[7];
	    aq = q2[7];
	    for (i__ = 6; i__ >= 0; --i__) {
		ap = p2[i__] + v * ap;
/* L2: */
		aq = q2[i__] + v * aq;
	    }
/* Computing 2nd power */
	    d__1 = v;
	    hc = exp(-(d__1 * d__1)) * ap / aq;
	    h__ = 1 - hc;
	} else {
/* Computing 2nd power */
	    d__1 = v;
	    y = 1 / (d__1 * d__1);
	    ap = p3[4];
	    aq = q3[4];
	    for (i__ = 3; i__ >= 0; --i__) {
		ap = p3[i__] + y * ap;
/* L3: */
		aq = q3[i__] + y * aq;
	    }
/* Computing 2nd power */
	    d__1 = v;
	    hc = exp(-(d__1 * d__1)) * (y * ap / aq + .56418958354775629) / v;
	    h__ = 1 - hc;
	}
	if (*x < 0.) {
	    h__ = -h__;
	    hc = 2 - hc;
	}
    }
    if (lef) {
	ret_val = h__;
    } else {
	ret_val = hc;
    }
    return ret_val;
} /* derf_ */

doublereal derf_(doublereal *x)
{
    return derf_0_(0, x);
    }

doublereal derfc_(doublereal *x)
{
    return derf_0_(1, x);
    }

