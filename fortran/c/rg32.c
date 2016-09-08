/* rg32.f -- translated by f2c (version 20090411).
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


/* $Id: rg32.F,v 1.1.1.1 1996/04/01 15:02:55 mclareni Exp $ */

/* $Log: rg32.F,v $ */
/* Revision 1.1.1.1  1996/04/01 15:02:55  mclareni */
/* Mathlib gen */


doublereal rg32_0_(int n__, integer *idummy)
{
    /* Initialized data */

    static integer iy = 875949887;

    /* System generated locals */
    real ret_val;

    /* Local variables */
    static integer jy, loop, jsum;

/*         MACHINE-INDEPENDENT RANDOM NUMBER GENERATOR */
/*         PRODUCES GAUSSIAN-DISTRIBUTED FLOATING-PT. */
/*               NUMBERS OF MEAN ZERO AND VARIANCE ONE. */
/*         IDENTICAL SEQUENCE ON ALL MACHINES OF .GE. 32 BITS. */
/*              UNIVERSAL VERSION,  F.JAMES, 1985 */
/*              IY     IS THE SEED, */
/*              CONS   IS 2**-23 */
/*              MASK31 IS 17777777777 OCTAL */
/* SELF,IF=-IBM. CERN EDITLIB+SIEMENS COMPILER PROBLEM. */
/* SELF. */
    switch(n__) {
	case 1: goto L_rg32in;
	case 2: goto L_rg32ot;
	}


    jsum = 0;
    for (loop = 1; loop <= 12; ++loop) {
	iy *= 69069;
/*         KEEP ONLY LOWER 31 BITS */
	iy &= 2147483647;
/*         SHIFT RIGHT 8 BITS TO AVOID OVFLO */
	jy = iy / 256;
	jsum += jy;
/* L50: */
    }
/*         CORRECT FOR TRUNCATION BIAS AND ROUND TO NEAREST */
/*         EVEN BYTE TO ASSURE EXACT FLOAT */
    jsum = (jsum + 134) / 256 << 8;

    ret_val = jsum * 1.1920928955078e-7f - 6.f;
    return ret_val;
/*         ENTRY TO INPUT SEED */

L_rg32in:
    iy = *idummy;
    ret_val = 0.f;
    return ret_val;
/*         ENTRY TO OUTPUT SEED */

L_rg32ot:
    *idummy = iy;
    ret_val = 0.f;
    return ret_val;
} /* rg32_ */

doublereal rg32_(integer *idummy)
{
    return rg32_0_(0, idummy);
    }

doublereal rg32in_(integer *idummy)
{
    return rg32_0_(1, idummy);
    }

doublereal rg32ot_(integer *idummy)
{
    return rg32_0_(2, idummy);
    }

