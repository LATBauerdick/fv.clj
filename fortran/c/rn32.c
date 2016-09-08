/* rn32.f -- translated by f2c (version 20090411).
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


/* $Id: rn32.F,v 1.1.1.1 1996/04/01 15:02:54 mclareni Exp $ */

/* $Log: rn32.F,v $ */
/* Revision 1.1.1.1  1996/04/01 15:02:54  mclareni */
/* Mathlib gen */


doublereal rn32_0_(int n__, integer *idummy)
{
    /* Initialized data */

    static integer iy = 65539;

    /* System generated locals */
    real ret_val;

    /* Local variables */
    static integer jy;

/*         MACHINE-INDEPENDENT RANDOM NUMBER GENERATOR */
/*         PRODUCES UNIFORMLY-DISTRIBUTED FLOATING-PT. */
/*               NUMBERS BETWEEN ZERO AND ONE. */
/*         IDENTICAL SEQUENCE ON ALL MACHINES OF .GE. 32 BITS. */
/*              UNIVERSAL VERSION,  F.JAMES, 1985 */
/*              IY     IS THE SEED, */
/*              CONS   IS 2**-31 */
/*              MASK31 IS 17777777777 OCTAL */
/* SELF,IF=-IBM. CERN EDITLIB+SIEMENS COMPILER PROBLEM. */
/* SELF. */
    switch(n__) {
	case 1: goto L_rn32in;
	case 2: goto L_rn32ot;
	}

    iy *= 69069;
/*         KEEP ONLY LOWER 31 BITS */
    iy &= 2147483647;
/*         SET LOWER 8 BITS TO ZERO TO ASSURE EXACT FLOAT */
    jy = iy / 256 << 8;
    ret_val = jy * 4.6566128730774e-10f;
    return ret_val;
/*         ENTRY TO INPUT SEED */

L_rn32in:
    iy = *idummy;
    ret_val = 0.f;
    return ret_val;
/*         ENTRY TO OUTPUT SEED */

L_rn32ot:
    *idummy = iy;
    ret_val = 0.f;
    return ret_val;
} /* rn32_ */

doublereal rn32_(integer *idummy)
{
    return rn32_0_(0, idummy);
    }

doublereal rn32in_(integer *idummy)
{
    return rn32_0_(1, idummy);
    }

doublereal rn32ot_(integer *idummy)
{
    return rn32_0_(2, idummy);
    }

