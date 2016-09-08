/* fvtmc.f -- translated by f2c (version 20090411).
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

/* Subroutine */ int fvtmcunify_(doublereal *p0, doublereal *p)
{
/* ----------------------------------------------------------------------- */
/* -- shift track parameters from p to p0 to a common vertex */
/* -- prelim.: by setting d0/z0 to zero */

    /* Parameter adjustments */
    --p;
    --p0;

    /* Function Body */
    p0[1] = p[1];
    p0[2] = p[2];
    p0[3] = p[3];
    p0[4] = .001;
    p0[5] = .001;
    return 0;
} /* fvtmcunify_ */

/* Subroutine */ int fvtmcsmear_(doublereal *p, doublereal *p0, doublereal *
	cplam, doublereal *cpu)
{
    static doublereal s;
    static integer i0, i1;
    static doublereal pp[5];
    extern doublereal rg32_(integer *);

/* ----------------------------------------------------------------------- */
/* -- smear track parameters p0 according to their error matrix */
/* -- denoted bu CpLam and CpU */

/* -- get vector in diagonalized system, pp = U^T.p0 */
    /* Parameter adjustments */
    cpu -= 6;
    --cplam;
    --p0;
    --p;

    /* Function Body */
    for (i0 = 1; i0 <= 5; ++i0) {
	s = 0.;
	for (i1 = 1; i1 <= 5; ++i1) {
	    s += cpu[i1 + i0 * 5] * p0[i1];
	}
	pp[i0 - 1] = s;
    }
/* -- smear pp with diag{CpLam} */
    for (i0 = 1; i0 <= 5; ++i0) {
	pp[i0 - 1] += cplam[i0] * (doublereal) rg32_(&i0);
    }
/* -- transform back p = U.pp */
    for (i0 = 1; i0 <= 5; ++i0) {
	s = 0.;
	for (i1 = 1; i1 <= 5; ++i1) {
	    s += cpu[i0 + i1 * 5] * pp[i1 - 1];
	}
	p[i0] = s;
    }
    return 0;
} /* fvtmcsmear_ */

