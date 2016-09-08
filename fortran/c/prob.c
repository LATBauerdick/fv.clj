/* prob.f -- translated by f2c (version 20090411).
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

/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b8 = .33333333333333331;


/* $Id: prob.F,v 1.1.1.1 1996/04/01 15:02:41 mclareni Exp $ */

/* $Log: prob.F,v $ */
/* Revision 1.1.1.1  1996/04/01 15:02:41  mclareni */
/* Mathlib gen */


doublereal prob_(real *x, integer *n)
{
    /* Format strings */
    static char fmt_101[] = "(\002PROB: N = \002,i6,\002 < 1\002)";
    static char fmt_102[] = "(\002PROB: X = \002,1p,e20.10,\002 < 0\002)";

    /* System generated locals */
    integer i__1;
    real ret_val;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;
    double sqrt(doublereal), pow_dd(doublereal *, doublereal *), exp(
	    doublereal);

    /* Local variables */
    static doublereal e, h__;
    static integer i__, m;
    static doublereal s, t, u, w, y, fi;
    extern doublereal derfc_(doublereal *);
    static char errtxt[80];
    extern /* Subroutine */ int fverror_(char *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___5 = { 0, errtxt, 0, fmt_101, 80, 1 };
    static icilist io___6 = { 0, errtxt, 0, fmt_102, 80, 1 };



/* $Id: imp64.inc,v 1.1.1.1 1996/04/01 15:02:59 mclareni Exp $ */

/* $Log: imp64.inc,v $ */
/* Revision 1.1.1.1  1996/04/01 15:02:59  mclareni */
/* Mathlib gen */


/* imp64.inc */

/*                maximum chi2 per df for df >= 2., if chi2/df > chipdf prob=0. */
    y = *x;
    u = y * .5;
    if (*n <= 0) {
	h__ = 0.;
	s_wsfi(&io___5);
	do_fio(&c__1, (char *)&(*n), (ftnlen)sizeof(integer));
	e_wsfi();
/*       CALL MTLPRT(NAME,'G100.1',ERRTXT) */
	fverror_(errtxt, (ftnlen)80);
    } else if (y < 0.) {
	h__ = 0.;
	s_wsfi(&io___6);
	do_fio(&c__1, (char *)&(*x), (ftnlen)sizeof(real));
	e_wsfi();
/*       CALL MTLPRT(NAME,'G100.2',ERRTXT) */
	fverror_(errtxt, (ftnlen)80);
    } else if (y == 0. || (doublereal) (*n / 20) > y) {
	h__ = 1.;
    } else if (*n == 1) {
	w = sqrt(u);
	if (w < 24.) {
	    h__ = derfc_(&w);
	} else {
	    h__ = 0.;
	}
    } else if (*n > 300) {
	s = 1. / *n;
	t = s * .22222222222222221;
	d__1 = y * s;
	w = (pow_dd(&d__1, &c_b8) - (1 - t)) / sqrt(t * 2);
	if (w < -24.) {
	    h__ = 1.;
	} else if (w < 24.) {
	    h__ = derfc_(&w) * .5;
	} else {
	    h__ = 0.;
	}
    } else {
	m = *n / 2;
	if (u < 349.346 && y / *n <= 100.) {
	    s = exp(u * -.5);
	    t = s;
	    e = s;
	    if (m << 1 == *n) {
		fi = 0.;
		i__1 = m - 1;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    fi += 1;
		    t = u * t / fi;
/* L1: */
		    s += t;
		}
		h__ = s * e;
	    } else {
		fi = 1.;
		i__1 = m - 1;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    fi += 2;
		    t = t * y / fi;
/* L2: */
		    s += t;
		}
		w = sqrt(u);
		if (w < 24.) {
		    h__ = w * 1.128379167095513 * s * e + derfc_(&w);
		} else {
		    h__ = 0.f;
		}
	    }
	} else {
	    h__ = 0.;
	}
    }
    if (h__ > 1e-30) {
	ret_val = h__;
    } else {
	ret_val = 0.f;
    }
    return ret_val;
} /* prob_ */

