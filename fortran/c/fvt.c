/* fvt.f -- translated by f2c (version 20090411).
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

/* Common Block Declarations */

struct {
    doublereal tw2pt, th[200]	/* was [5][40] */, tch[1000]	/* was [5][5][
	    40] */, tx[120]	/* was [3][40] */, tcx[360]	/* was [3][3][
	    40] */, tgh[1000]	/* was [5][5][40] */, tt[160]	/* was [4][40]
	     */, tct[640]	/* was [4][4][40] */, tq[120]	/* was [3][40]
	     */, tcq[360]	/* was [3][3][40] */, tchi2[40];
    integer tnt;
    logical fvtprint;
    integer fvtplun;
} fvtcom_;

#define fvtcom_1 fvtcom_

/* Table of constant values */

static integer c__1 = 1;
static integer c__6 = 6;
static integer c__5 = 5;
static integer c__4 = 4;
static integer c__2 = 2;
static integer c__3 = 3;
static integer c__9 = 9;

/* ---------------------------------------------------------------------- */

/* FVT, test the FV package to fit vertices and find secondary vertices */

/* author: Lothar A.T. Bauerdick, CERN/PPE */
/* date:   Jan. 1991 */

/* ref.:   R. Fruehwirt, ``Applications of Filter Methods...'' */
/*         HEPHY-PUB 516/88, Vienna */

/* ---------------------------------------------------------------------- */
/* Main program */ int MAIN__(void)
{
    /* System generated locals */
    olist o__1;

    /* Builtin functions */
    integer f_open(olist *);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer do__;
    extern /* Subroutine */ int dom5_(void), dowu_(void), doors_(void), 
	    dom5mc_(void), fvsprt_(logical *, integer *), errinit_(void), 
	    errsumm_(void);

/* CC      parameter (dataDir = 'disk$user:[bauerdick.fv.dat]') */
    do__ = 1;
    fvtcom_1.fvtprint = TRUE_;
    fvtcom_1.fvtplun = 41;
    fvsprt_(&fvtcom_1.fvtprint, &fvtcom_1.fvtplun);
/* -- default value for mag.field */
    fvtcom_1.tw2pt = .0044968501;
    errinit_();
    if (fvtcom_1.fvtprint) {
	o__1.oerr = 0;
	o__1.ounit = fvtcom_1.fvtplun;
	o__1.ofnmlen = 10;
	o__1.ofnm = "fvt.output";
	o__1.orl = 0;
	o__1.osta = "UNKNOWN";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
    }
    if (do__ == 5) {
/* CCCCCCCCCCCCC        call doMc */
    } else if (do__ == 1) {
	dom5_();
    } else if (do__ == 2) {
	dom5mc_();
    } else if (do__ == 3) {
/* CCCCCCCCCCCCC        call doHad */
    } else if (do__ == 4) {
	dowu_();
    } else {
	doors_();
    }
    errsumm_();
    s_stop("", (ftnlen)0);
    return 0;
} /* MAIN__ */

/* Subroutine */ int dom5mc_(void)
{
    /* Initialized data */

    static integer runlist[6]	/* was [2][3] */ = { 101,7076,101,8340,101,
	    12040 };
    static integer trlist[24]	/* was [6][4] */ = { 1,2,3,4,5,6,1,3,4,5,6,2,
	    1,2,3,5,6,4,1,3,4,5,6,2 };
    static integer plist[24]	/* was [6][4] */ = { 4,4,4,4,4,4,4,4,4,4,4,4,
	    4,4,4,4,4,4,4,4,4,4,4,4 };

    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static integer i__;
    static doublereal m;
    static integer i0, i1;
    extern /* Subroutine */ int fvthelix2p4_(doublereal *, doublereal *, 
	    integer *, integer *, integer *, doublereal *, doublereal *);
    static doublereal chi2, sigm;
    static integer kevt, krun;
    extern integer fvtdo_(doublereal *, integer *, integer *);
    static integer status;
    extern integer fvtread_(integer *, integer *);
    extern /* Subroutine */ int fvtqm2p4_(doublereal *, doublereal *, integer 
	    *, integer *, integer *, doublereal *, doublereal *);
    extern integer fvtinvm_(doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *);

    /* Fortran I/O blocks */
    static cilist io___8 = { 0, 0, 0, "(1x,78(1H=))", 0 };
    static cilist io___9 = { 0, 0, 0, "(/1x,a,i5,a,i5)", 0 };
    static cilist io___14 = { 0, 0, 0, "(/1x,a,i5,a,i5)", 0 };
    static cilist io___17 = { 0, 0, 0, "(1x,a,g12.5,a,g12.5)", 0 };
    static cilist io___18 = { 0, 0, 0, "(1x,a,g12.5,a,g12.5)", 0 };


/* -------------------------------------------------------------------- */
/* CC      parameter (dataDir = 'disk$user:[bauerdick.fv.dat]') */
/* -- 1:electron, 2:muon, 3:tau, 4:pi */
    for (i__ = 1; i__ <= 3; ++i__) {
	krun = runlist[(i__ << 1) - 2];
	kevt = runlist[(i__ << 1) - 1];
	if (fvtcom_1.fvtprint) {
	    io___8.ciunit = fvtcom_1.fvtplun;
	    s_wsfe(&io___8);
	    e_wsfe();
	    io___9.ciunit = fvtcom_1.fvtplun;
	    s_wsfe(&io___9);
	    do_fio(&c__1, "Run ", (ftnlen)4);
	    do_fio(&c__1, (char *)&krun, (ftnlen)sizeof(integer));
	    do_fio(&c__1, " Evt ", (ftnlen)5);
	    do_fio(&c__1, (char *)&kevt, (ftnlen)sizeof(integer));
	    e_wsfe();
	}
	status = fvtread_(&krun, &kevt);
	for (i0 = 1; i0 <= 3; ++i0) {
	    for (i1 = 1; i1 <= 3; ++i1) {
		fvtcom_1.tcx[i0 + (i1 + 3) * 3 - 13] *= 1e4;
	    }
	}
	status = fvtdo_(&chi2, &c__6, &trlist[i__ * 6]);
/* -- */
	status = fvtread_(&krun, &kevt);
	for (i0 = 1; i0 <= 3; ++i0) {
	    for (i1 = 1; i1 <= 3; ++i1) {
		fvtcom_1.tcx[i0 + (i1 + 3) * 3 - 13] *= 1e4;
	    }
	}
	if (fvtcom_1.fvtprint) {
	    io___14.ciunit = fvtcom_1.fvtplun;
	    s_wsfe(&io___14);
	    do_fio(&c__1, "common vertex for 5 tracks for Run ", (ftnlen)35);
	    do_fio(&c__1, (char *)&krun, (ftnlen)sizeof(integer));
	    do_fio(&c__1, " Evt ", (ftnlen)5);
	    do_fio(&c__1, (char *)&kevt, (ftnlen)sizeof(integer));
	    e_wsfe();
/* -- fill tt, tCt from th, tCh for particles trList, pList */
/* -- and calculate inv. mass from  non-fitted momenta */
	    fvthelix2p4_(fvtcom_1.tt, fvtcom_1.tct, &c__5, &trlist[i__ * 6], &
		    plist[i__ * 6], fvtcom_1.th, fvtcom_1.tch);
	    status = fvtinvm_(&m, &sigm, &c__5, &trlist[i__ * 6], fvtcom_1.tt,
		     fvtcom_1.tct);
	    io___17.ciunit = fvtcom_1.fvtplun;
	    s_wsfe(&io___17);
	    do_fio(&c__1, "invariant mass is ", (ftnlen)18);
	    d__1 = m * 1000;
	    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, " +/- ", (ftnlen)5);
	    d__2 = sigm * 1000;
	    do_fio(&c__1, (char *)&d__2, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
	status = fvtdo_(&chi2, &c__5, &trlist[i__ * 6]);
	if (fvtcom_1.fvtprint) {
/* -- calculate Energy for particle types in pList */
	    fvtqm2p4_(fvtcom_1.tt, fvtcom_1.tct, &c__5, &trlist[i__ * 6], &
		    plist[i__ * 6], fvtcom_1.tq, fvtcom_1.tcq);
	    status = fvtinvm_(&m, &sigm, &c__5, &trlist[i__ * 6], fvtcom_1.tt,
		     fvtcom_1.tct);
	    io___18.ciunit = fvtcom_1.fvtplun;
	    s_wsfe(&io___18);
	    do_fio(&c__1, "invariant mass is ", (ftnlen)18);
	    d__1 = m * 1000;
	    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, " +/- ", (ftnlen)5);
	    d__2 = sigm * 1000;
	    do_fio(&c__1, (char *)&d__2, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }
    return 0;
} /* dom5mc_ */

/* Subroutine */ int dom5_(void)
{
    /* Initialized data */

    static integer runlist[16]	/* was [2][8] */ = { 5129,1412,5158,4656,5166,
	    1984,5343,2291,7849,7984,8489,4451,8985,5486,9019,4769 };
    static integer trlist[54]	/* was [6][9] */ = { 1,2,3,4,5,6,1,3,4,5,6,2,
	    1,2,3,5,6,4,2,3,4,5,6,1,1,3,4,5,6,2,1,3,4,5,6,2,1,2,3,5,6,4,1,2,3,
	    4,6,5,1,2,3,4,6,5 };
    static integer plist[54]	/* was [6][9] */ = { 4,4,4,4,4,4,4,4,4,4,4,4,
	    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
	    4,4,4,4,4,4,4,4,4 };

    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static integer i__;
    static doublereal m;
    static integer i0, i1;
    extern /* Subroutine */ int fvthelix2p4_(doublereal *, doublereal *, 
	    integer *, integer *, integer *, doublereal *, doublereal *);
    static doublereal chi2, sigm;
    static integer kevt, krun;
    extern integer fvtdo_(doublereal *, integer *, integer *);
    static integer status;
    extern integer fvtread_(integer *, integer *);
    extern /* Subroutine */ int fvtqm2p4_(doublereal *, doublereal *, integer 
	    *, integer *, integer *, doublereal *, doublereal *);
    extern integer fvtinvm_(doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *);

    /* Fortran I/O blocks */
    static cilist io___25 = { 0, 0, 0, "(1x,78(1H=))", 0 };
    static cilist io___26 = { 0, 0, 0, "(/1x,a,i5,a,i5)", 0 };
    static cilist io___31 = { 0, 0, 0, "(/1x,a,i5,a,i5)", 0 };
    static cilist io___34 = { 0, 0, 0, "(1x,a,g12.5,a,g12.5)", 0 };
    static cilist io___35 = { 0, 0, 0, "(1x,a,g12.5,a,g12.5)", 0 };


/* -------------------------------------------------------------------- */
/* CC      parameter (dataDir = 'disk$user:[bauerdick.fv.dat]') */
/* -- 1:electron, 2:muon, 3:tau, 4:pi */
    for (i__ = 1; i__ <= 8; ++i__) {
	krun = runlist[(i__ << 1) - 2];
	kevt = runlist[(i__ << 1) - 1];
	if (fvtcom_1.fvtprint) {
	    io___25.ciunit = fvtcom_1.fvtplun;
	    s_wsfe(&io___25);
	    e_wsfe();
	    io___26.ciunit = fvtcom_1.fvtplun;
	    s_wsfe(&io___26);
	    do_fio(&c__1, "Run ", (ftnlen)4);
	    do_fio(&c__1, (char *)&krun, (ftnlen)sizeof(integer));
	    do_fio(&c__1, " Evt ", (ftnlen)5);
	    do_fio(&c__1, (char *)&kevt, (ftnlen)sizeof(integer));
	    e_wsfe();
	}
	status = fvtread_(&krun, &kevt);
	for (i0 = 1; i0 <= 3; ++i0) {
	    for (i1 = 1; i1 <= 3; ++i1) {
		fvtcom_1.tcx[i0 + (i1 + 3) * 3 - 13] *= 1e4;
	    }
	}
	status = fvtdo_(&chi2, &c__6, &trlist[i__ * 6]);
	status = fvtread_(&krun, &kevt);
	for (i0 = 1; i0 <= 3; ++i0) {
	    for (i1 = 1; i1 <= 3; ++i1) {
		fvtcom_1.tcx[i0 + (i1 + 3) * 3 - 13] *= 1e4;
	    }
	}
	if (fvtcom_1.fvtprint) {
	    io___31.ciunit = fvtcom_1.fvtplun;
	    s_wsfe(&io___31);
	    do_fio(&c__1, "common vertex for 5 tracks for Run ", (ftnlen)35);
	    do_fio(&c__1, (char *)&krun, (ftnlen)sizeof(integer));
	    do_fio(&c__1, " Evt ", (ftnlen)5);
	    do_fio(&c__1, (char *)&kevt, (ftnlen)sizeof(integer));
	    e_wsfe();
/* -- fill tt, tCt from th, tCh for particles trList, pList */
	    fvthelix2p4_(fvtcom_1.tt, fvtcom_1.tct, &c__5, &trlist[i__ * 6], &
		    plist[i__ * 6], fvtcom_1.th, fvtcom_1.tch);
/* -- calculate inv. mass */
	    status = fvtinvm_(&m, &sigm, &c__5, &trlist[i__ * 6], fvtcom_1.tt,
		     fvtcom_1.tct);
	    io___34.ciunit = fvtcom_1.fvtplun;
	    s_wsfe(&io___34);
	    do_fio(&c__1, "invariant mass is ", (ftnlen)18);
	    d__1 = m * 1000;
	    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, " +/- ", (ftnlen)5);
	    d__2 = sigm * 1000;
	    do_fio(&c__1, (char *)&d__2, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
	status = fvtdo_(&chi2, &c__5, &trlist[i__ * 6]);
	if (fvtcom_1.fvtprint) {
/* -- calculate Energy for particle types in pList */
	    fvtqm2p4_(fvtcom_1.tt, fvtcom_1.tct, &c__5, &trlist[i__ * 6], &
		    plist[i__ * 6], fvtcom_1.tq, fvtcom_1.tcq);
	    status = fvtinvm_(&m, &sigm, &c__5, &trlist[i__ * 6], fvtcom_1.tt,
		     fvtcom_1.tct);
	    io___35.ciunit = fvtcom_1.fvtplun;
	    s_wsfe(&io___35);
	    do_fio(&c__1, "invariant mass is ", (ftnlen)18);
	    d__1 = m * 1000;
	    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, " +/- ", (ftnlen)5);
	    d__2 = sigm * 1000;
	    do_fio(&c__1, (char *)&d__2, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }
    return 0;
} /* dom5_ */

/* Subroutine */ int dowu_(void)
{
    /* Initialized data */

    static integer runlist[32]	/* was [2][16] */ = { 5166,981,5825,51,6972,
	    4547,7422,2873,7496,3332,7681,2160,8126,771,8163,2748,8334,7367,
	    8335,1133,8537,4398,8623,158,8833,4673,8855,3264,8865,6111,8960,
	    749 };
    static integer trlist[68]	/* was [4][17] */ = { 1,2,3,4,4,3,2,1,3,4,1,2,
	    1,4,2,3,3,2,4,1,1,2,4,3,4,1,2,3,3,4,1,2,3,2,1,4,2,4,1,3,2,3,4,1,3,
	    4,2,1,4,2,3,1,3,2,4,1,4,1,3,2,1,2,4,3,4,1,2,3 };
    static integer plist[68]	/* was [4][17] */ = { 4,4,4,4,2,2,4,2,1,1,2,2,
	    4,4,4,4,4,4,1,1,4,4,2,2,1,1,1,4,4,4,1,4,4,4,1,1,4,4,1,1,1,1,4,1,4,
	    4,2,4,2,2,4,4,4,4,4,4,1,1,2,2,2,2,2,1,1,1,1,1 };
    static doublereal mass[4] = { 5.1099906e-4,.105658387,1.7841,.1395675 };

    /* System generated locals */
    real r__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen),
	     s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);

    /* Local variables */
    static doublereal d__;
    static integer i__;
    static doublereal m, p[3], x[3];
    static integer i0, i1;
    static doublereal dp[9]	/* was [3][3] */, ep[9]	/* was [3][3] */, ml[
	    5], cx[9]	/* was [3][3] */, pr;
    extern /* Subroutine */ int fvthelix2p4_(doublereal *, doublereal *, 
	    integer *, integer *, integer *, doublereal *, doublereal *);
    static integer new__;
    static doublereal chi2, sigd, sigm;
    extern doublereal prob_(real *, integer *);
    static integer kevt, krun;
    extern integer fvtdo_(doublereal *, integer *, integer *), fvfitd_(
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    extern /* Subroutine */ int fvdist_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), fvinvm_(
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *), fvcopy_(doublereal *, doublereal *, 
	    integer *, integer *);
    extern integer fvsumq_(doublereal *, doublereal *, doublereal *, integer *
	    , integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    static integer status;
    extern integer fvtread_(integer *, integer *);
    extern /* Subroutine */ int fvtqm2p4_(doublereal *, doublereal *, integer 
	    *, integer *, integer *, doublereal *, doublereal *);
    extern integer fvtinvm_(doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *);

    /* Fortran I/O blocks */
    static cilist io___43 = { 0, 0, 0, "(1x,78(1H=))", 0 };
    static cilist io___44 = { 0, 0, 0, "(/1x,a,i5,a,i5)", 0 };
    static cilist io___49 = { 0, 0, 0, "(/1x,a,i5,a,i5)", 0 };
    static cilist io___52 = { 0, 0, 0, "(1x,a,g12.5,a,g12.5)", 0 };
    static cilist io___56 = { 0, 0, 0, "(1x,a,g12.5,a,g12.5)", 0 };
    static cilist io___61 = { 0, 0, 0, "(/1x,a,i5,a,i5)", 0 };
    static cilist io___62 = { 0, 0, 0, "(1x,a,g12.5,a,g12.5)", 0 };
    static cilist io___63 = { 0, 0, 0, "(1x,a,g12.5,a,g12.5)", 0 };
    static cilist io___64 = { 0, 0, 0, "(1x,a,g12.5,a,g12.5)", 0 };
    static cilist io___67 = { 0, 0, 0, "(/1x,a,g10.3,a,g10.3)", 0 };
    static cilist io___68 = { 0, 6, 0, 0, 0 };
    static cilist io___70 = { 0, 0, 0, "(1x,a,g10.3,a,g10.3,a,g10.3,a,g10.3)",
	     0 };


/* -------------------------------------------------------------------- */
/* CC      parameter (dataDir = 'disk$user:[bauerdick.fv.dat]') */
/* -- tracks: */
/* -- event: Z0 -> l0 l2 */
/* --              l0 -> l1 V+ V- */
/* -- V+, V-, l1, l2 */
/* -- 1:electron, 2:muon, 3:tau, 4:pi */
    for (i__ = 1; i__ <= 16; ++i__) {
	krun = runlist[(i__ << 1) - 2];
	kevt = runlist[(i__ << 1) - 1];
/* -- fit for common vertex of l+l-V=V- */
	io___43.ciunit = fvtcom_1.fvtplun;
	s_wsfe(&io___43);
	e_wsfe();
	io___44.ciunit = fvtcom_1.fvtplun;
	s_wsfe(&io___44);
	do_fio(&c__1, "Run ", (ftnlen)4);
	do_fio(&c__1, (char *)&krun, (ftnlen)sizeof(integer));
	do_fio(&c__1, " Evt ", (ftnlen)5);
	do_fio(&c__1, (char *)&kevt, (ftnlen)sizeof(integer));
	e_wsfe();
	status = fvtread_(&krun, &kevt);
	for (i0 = 1; i0 <= 3; ++i0) {
	    for (i1 = 1; i1 <= 3; ++i1) {
		fvtcom_1.tcx[i0 + (i1 + 3) * 3 - 13] *= 1e4;
	    }
	}
	status = fvtdo_(&chi2, &c__4, trlist);
/* -- fit for vertex of V+V- */
	status = fvtread_(&krun, &kevt);
	for (i0 = 1; i0 <= 3; ++i0) {
	    for (i1 = 1; i1 <= 3; ++i1) {
		fvtcom_1.tcx[i0 + (i1 + 3) * 3 - 13] *= 1e4;
	    }
	}
	io___49.ciunit = fvtcom_1.fvtplun;
	s_wsfe(&io___49);
	do_fio(&c__1, "common vertex for V+V- for Run ", (ftnlen)31);
	do_fio(&c__1, (char *)&krun, (ftnlen)sizeof(integer));
	do_fio(&c__1, " Evt ", (ftnlen)5);
	do_fio(&c__1, (char *)&kevt, (ftnlen)sizeof(integer));
	e_wsfe();
/* -- fill tt, tCt from th, tCh for particles trList, pList */
	fvthelix2p4_(fvtcom_1.tt, fvtcom_1.tct, &c__2, &trlist[i__ * 4], &
		plist[i__ * 4], fvtcom_1.th, fvtcom_1.tch);
	status = fvtinvm_(&m, &sigm, &c__2, &trlist[i__ * 4], fvtcom_1.tt, 
		fvtcom_1.tct);
	io___52.ciunit = fvtcom_1.fvtplun;
	s_wsfe(&io___52);
	do_fio(&c__1, "invariant mass is ", (ftnlen)18);
	d__1 = m * 1000;
	do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, " +/- ", (ftnlen)5);
	d__2 = sigm * 1000;
	do_fio(&c__1, (char *)&d__2, (ftnlen)sizeof(doublereal));
	e_wsfe();
	status = fvtdo_(&chi2, &c__2, &trlist[i__ * 4]);
/* -- save the V+V- vertex in x, Cx */
	fvcopy_(x, fvtcom_1.tx, &c__3, &c__1);
	fvcopy_(cx, fvtcom_1.tcx, &c__3, &c__3);
/* -- calculate invariant mass of V+V- */
	ml[trlist[i__ * 4] - 1] = mass[plist[i__ * 4] - 1];
	ml[trlist[(i__ << 2) + 1] - 1] = mass[plist[(i__ << 2) + 1] - 1];
	fvinvm_(&m, &sigm, &c__2, &trlist[i__ * 4], &fvtcom_1.tw2pt, ml, x, 
		cx, fvtcom_1.tq, fvtcom_1.tcq, fvtcom_1.tgh);
	io___56.ciunit = fvtcom_1.fvtplun;
	s_wsfe(&io___56);
	do_fio(&c__1, "invariant mass is ", (ftnlen)18);
	d__1 = m * 1000;
	do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, " +/- ", (ftnlen)5);
	d__2 = sigm * 1000;
	do_fio(&c__1, (char *)&d__2, (ftnlen)sizeof(doublereal));
	e_wsfe();
/* -- sum up momenta of V+V- and put resultant track parameters into th/tGh */
	new__ = fvtcom_1.tnt + 1;
	status = fvsumq_(p, dp, ep, &new__, &fvtcom_1.tnt, trlist, &
		fvtcom_1.tw2pt, x, cx, fvtcom_1.tq, fvtcom_1.tcq, fvtcom_1.th,
		 fvtcom_1.tgh);
/* -- fit vertex to l+l- pair */
	io___61.ciunit = fvtcom_1.fvtplun;
	s_wsfe(&io___61);
	do_fio(&c__1, "common vertex for l(V)l for Run ", (ftnlen)32);
	do_fio(&c__1, (char *)&krun, (ftnlen)sizeof(integer));
	do_fio(&c__1, " Evt ", (ftnlen)5);
	do_fio(&c__1, (char *)&kevt, (ftnlen)sizeof(integer));
	e_wsfe();
	fvthelix2p4_(fvtcom_1.tt, fvtcom_1.tct, &c__2, &trlist[(i__ << 2) + 2]
		, &plist[(i__ << 2) + 2], fvtcom_1.th, fvtcom_1.tch);
	status = fvtinvm_(&m, &sigm, &c__2, &trlist[(i__ << 2) + 2], 
		fvtcom_1.tt, fvtcom_1.tct);
	io___62.ciunit = fvtcom_1.fvtplun;
	s_wsfe(&io___62);
	do_fio(&c__1, "invariant mass is ", (ftnlen)18);
	d__1 = m * 1000;
	do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, " +/- ", (ftnlen)5);
	d__2 = sigm * 1000;
	do_fio(&c__1, (char *)&d__2, (ftnlen)sizeof(doublereal));
	e_wsfe();
	for (i0 = 1; i0 <= 3; ++i0) {
	    for (i1 = 1; i1 <= 3; ++i1) {
		fvtcom_1.tcx[i0 + (i1 + 3) * 3 - 13] *= 1e4;
	    }
	}
	status = fvtdo_(&chi2, &c__2, &trlist[(i__ << 2) + 2]);
	ml[trlist[(i__ << 2) + 2] - 1] = mass[plist[(i__ << 2) + 2] - 1];
	ml[trlist[(i__ << 2) + 3] - 1] = mass[plist[(i__ << 2) + 3] - 1];
	fvinvm_(&m, &sigm, &c__2, &trlist[(i__ << 2) + 2], &fvtcom_1.tw2pt, 
		ml, fvtcom_1.tx, fvtcom_1.tcx, fvtcom_1.tq, fvtcom_1.tcq, 
		fvtcom_1.tgh);
	io___63.ciunit = fvtcom_1.fvtplun;
	s_wsfe(&io___63);
	do_fio(&c__1, "invariant mass is ", (ftnlen)18);
	d__1 = m * 1000;
	do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, " +/- ", (ftnlen)5);
	d__2 = sigm * 1000;
	do_fio(&c__1, (char *)&d__2, (ftnlen)sizeof(doublereal));
	e_wsfe();
	fvtqm2p4_(fvtcom_1.tt, fvtcom_1.tct, &c__2, &trlist[(i__ << 2) + 2], &
		plist[(i__ << 2) + 2], fvtcom_1.tq, fvtcom_1.tcq);
	status = fvtinvm_(&m, &sigm, &c__2, &trlist[(i__ << 2) + 2], 
		fvtcom_1.tt, fvtcom_1.tct);
	io___64.ciunit = fvtcom_1.fvtplun;
	s_wsfe(&io___64);
	do_fio(&c__1, "invariant mass is ", (ftnlen)18);
	d__1 = m * 1000;
	do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, " +/- ", (ftnlen)5);
	d__2 = sigm * 1000;
	do_fio(&c__1, (char *)&d__2, (ftnlen)sizeof(doublereal));
	e_wsfe();
/* -- calculate distance between vertex of V+V- and l+l- */
	fvdist_(&d__, &sigd, x, cx, fvtcom_1.tx, fvtcom_1.tcx);
	io___67.ciunit = fvtcom_1.fvtplun;
	s_wsfe(&io___67);
	do_fio(&c__1, "distance between V/ll vertices is ", (ftnlen)34);
	do_fio(&c__1, (char *)&d__, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, " +/- ", (ftnlen)5);
	do_fio(&c__1, (char *)&sigd, (ftnlen)sizeof(doublereal));
	e_wsfe();
/* -- fit distance between the two */
	status = fvfitd_(&d__, &sigd, &chi2, &d__, fvtcom_1.tx, fvtcom_1.tcx, 
		x, cx, p, dp, ep);
	if ((status & 1) != 1) {
	    s_wsle(&io___68);
	    do_lio(&c__9, &c__1, "fit to tau did not work...", (ftnlen)26);
	    e_wsle();
	}
	r__1 = (real) chi2;
	pr = prob_(&r__1, &c__1);
	io___70.ciunit = fvtcom_1.fvtplun;
	s_wsfe(&io___70);
	do_fio(&c__1, "fitted distance between V/ll vertices is ", (ftnlen)41)
		;
	do_fio(&c__1, (char *)&d__, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, " +/- ", (ftnlen)5);
	do_fio(&c__1, (char *)&sigd, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, " chi^2 ", (ftnlen)7);
	do_fio(&c__1, (char *)&chi2, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, " prob. ", (ftnlen)7);
	do_fio(&c__1, (char *)&pr, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    return 0;
} /* dowu_ */

/* Subroutine */ int doors_(void)
{
    /* Initialized data */

    static integer runlist[22]	/* was [2][11] */ = { 4515,750,5060,2137,5821,
	    3675,5880,3296,7339,1745,7412,4101,7743,8370,8335,1133,8383,966,
	    8619,2198,8898,5940 };
    static integer trlist[48]	/* was [4][12] */ = { 1,2,3,4,3,4,1,2,4,1,3,2,
	    2,4,3,1,1,4,2,3,4,3,2,1,4,3,2,1,5,6,1,3,4,3,2,1,5,6,1,2,3,4,2,1,3,
	    4,2,5 };
    static integer plist[48]	/* was [4][12] */ = { 4,4,4,4,4,4,2,2,4,4,1,1,
	    2,2,2,2,1,1,1,1,1,1,2,2,1,1,2,2,1,1,4,4,1,1,4,4,1,1,4,4,1,1,2,2,1,
	    1,1,1 };

    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static doublereal d__;
    static integer i__;
    static doublereal m, x[3];
    static integer i0, i1;
    static doublereal cx[9]	/* was [3][3] */;
    extern /* Subroutine */ int fvthelix2p4_(doublereal *, doublereal *, 
	    integer *, integer *, integer *, doublereal *, doublereal *);
    static doublereal chi2, sigd, sigm;
    static integer kevt, krun;
    extern integer fvtdo_(doublereal *, integer *, integer *);
    extern /* Subroutine */ int fvdist_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), fvcopy_(
	    doublereal *, doublereal *, integer *, integer *);
    static integer status;
    extern integer fvtread_(integer *, integer *);
    extern /* Subroutine */ int fvtqm2p4_(doublereal *, doublereal *, integer 
	    *, integer *, integer *, doublereal *, doublereal *);
    extern integer fvtinvm_(doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *);

    /* Fortran I/O blocks */
    static cilist io___77 = { 0, 0, 0, "(1x,78(1H=))", 0 };
    static cilist io___78 = { 0, 0, 0, "(/1x,a,i5,a,i5)", 0 };
    static cilist io___83 = { 0, 0, 0, "(/1x,a,i5,a,i5)", 0 };
    static cilist io___86 = { 0, 0, 0, "(1x,a,g12.5,a,g12.5)", 0 };
    static cilist io___89 = { 0, 0, 0, "(1x,a,g12.5,a,g12.5)", 0 };
    static cilist io___90 = { 0, 0, 0, "(/1x,a,i5,a,i5)", 0 };
    static cilist io___91 = { 0, 0, 0, "(1x,a,g12.5,a,g12.5)", 0 };
    static cilist io___92 = { 0, 0, 0, "(1x,a,g12.5,a,g12.5)", 0 };
    static cilist io___95 = { 0, 0, 0, "(/1x,a,g10.3,a,g10.3)", 0 };


/* -------------------------------------------------------------------- */
/* CC      parameter (dataDir = 'disk$user:[bauerdick.fv.dat]') */
/* -- tracks: */
/* -- event: Z0 -> l0 l2 */
/* --              l0 -> l1 V+ V- */
/* -- V+, V-, l1, l2 */
/* -- 1:electron, 2:muon, 3:tau, 4:pi */
    for (i__ = 1; i__ <= 11; ++i__) {
	krun = runlist[(i__ << 1) - 2];
	kevt = runlist[(i__ << 1) - 1];
	io___77.ciunit = fvtcom_1.fvtplun;
	s_wsfe(&io___77);
	e_wsfe();
	io___78.ciunit = fvtcom_1.fvtplun;
	s_wsfe(&io___78);
	do_fio(&c__1, "Run ", (ftnlen)4);
	do_fio(&c__1, (char *)&krun, (ftnlen)sizeof(integer));
	do_fio(&c__1, " Evt ", (ftnlen)5);
	do_fio(&c__1, (char *)&kevt, (ftnlen)sizeof(integer));
	e_wsfe();
	status = fvtread_(&krun, &kevt);
	for (i0 = 1; i0 <= 3; ++i0) {
	    for (i1 = 1; i1 <= 3; ++i1) {
		fvtcom_1.tcx[i0 + (i1 + 3) * 3 - 13] *= 1e4;
	    }
	}
	status = fvtdo_(&chi2, &c__4, trlist);
	status = fvtread_(&krun, &kevt);
	for (i0 = 1; i0 <= 3; ++i0) {
	    for (i1 = 1; i1 <= 3; ++i1) {
		fvtcom_1.tcx[i0 + (i1 + 3) * 3 - 13] *= 1e4;
	    }
	}
	io___83.ciunit = fvtcom_1.fvtplun;
	s_wsfe(&io___83);
	do_fio(&c__1, "common vertex for V+V- for Run ", (ftnlen)31);
	do_fio(&c__1, (char *)&krun, (ftnlen)sizeof(integer));
	do_fio(&c__1, " Evt ", (ftnlen)5);
	do_fio(&c__1, (char *)&kevt, (ftnlen)sizeof(integer));
	e_wsfe();
/* -- fill tt, tCt from th, tCh for particles trList, pList */
	fvthelix2p4_(fvtcom_1.tt, fvtcom_1.tct, &c__2, &trlist[i__ * 4], &
		plist[i__ * 4], fvtcom_1.th, fvtcom_1.tch);
	status = fvtinvm_(&m, &sigm, &c__2, &trlist[i__ * 4], fvtcom_1.tt, 
		fvtcom_1.tct);
	io___86.ciunit = fvtcom_1.fvtplun;
	s_wsfe(&io___86);
	do_fio(&c__1, "invariant mass is ", (ftnlen)18);
	d__1 = m * 1000;
	do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, " +/- ", (ftnlen)5);
	d__2 = sigm * 1000;
	do_fio(&c__1, (char *)&d__2, (ftnlen)sizeof(doublereal));
	e_wsfe();
	status = fvtdo_(&chi2, &c__2, &trlist[i__ * 4]);
	fvcopy_(x, fvtcom_1.tx, &c__3, &c__1);
	fvcopy_(cx, fvtcom_1.tcx, &c__3, &c__3);
	fvtqm2p4_(fvtcom_1.tt, fvtcom_1.tct, &c__2, &trlist[i__ * 4], &plist[
		i__ * 4], fvtcom_1.tq, fvtcom_1.tcq);
	status = fvtinvm_(&m, &sigm, &c__2, &trlist[i__ * 4], fvtcom_1.tt, 
		fvtcom_1.tct);
	io___89.ciunit = fvtcom_1.fvtplun;
	s_wsfe(&io___89);
	do_fio(&c__1, "invariant mass is ", (ftnlen)18);
	d__1 = m * 1000;
	do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, " +/- ", (ftnlen)5);
	d__2 = sigm * 1000;
	do_fio(&c__1, (char *)&d__2, (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___90.ciunit = fvtcom_1.fvtplun;
	s_wsfe(&io___90);
	do_fio(&c__1, "common vertex for l(V)l for Run ", (ftnlen)32);
	do_fio(&c__1, (char *)&krun, (ftnlen)sizeof(integer));
	do_fio(&c__1, " Evt ", (ftnlen)5);
	do_fio(&c__1, (char *)&kevt, (ftnlen)sizeof(integer));
	e_wsfe();
	fvthelix2p4_(fvtcom_1.tt, fvtcom_1.tct, &c__2, &trlist[(i__ << 2) + 2]
		, &plist[i__ * 4], fvtcom_1.th, fvtcom_1.tch);
	status = fvtinvm_(&m, &sigm, &c__2, &trlist[(i__ << 2) + 2], 
		fvtcom_1.tt, fvtcom_1.tct);
	io___91.ciunit = fvtcom_1.fvtplun;
	s_wsfe(&io___91);
	do_fio(&c__1, "invariant mass is ", (ftnlen)18);
	d__1 = m * 1000;
	do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, " +/- ", (ftnlen)5);
	d__2 = sigm * 1000;
	do_fio(&c__1, (char *)&d__2, (ftnlen)sizeof(doublereal));
	e_wsfe();
	status = fvtdo_(&chi2, &c__2, &trlist[(i__ << 2) + 2]);
	fvtqm2p4_(fvtcom_1.tt, fvtcom_1.tct, &c__2, &trlist[(i__ << 2) + 2], &
		plist[i__ * 4], fvtcom_1.tq, fvtcom_1.tcq);
	status = fvtinvm_(&m, &sigm, &c__2, &trlist[(i__ << 2) + 2], 
		fvtcom_1.tt, fvtcom_1.tct);
	io___92.ciunit = fvtcom_1.fvtplun;
	s_wsfe(&io___92);
	do_fio(&c__1, "invariant mass is ", (ftnlen)18);
	d__1 = m * 1000;
	do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, " +/- ", (ftnlen)5);
	d__2 = sigm * 1000;
	do_fio(&c__1, (char *)&d__2, (ftnlen)sizeof(doublereal));
	e_wsfe();
/* -- calculate distance */
	fvdist_(&d__, &sigd, x, cx, fvtcom_1.tx, fvtcom_1.tcx);
	io___95.ciunit = fvtcom_1.fvtplun;
	s_wsfe(&io___95);
	do_fio(&c__1, "distance between V/ll vertices is ", (ftnlen)34);
	do_fio(&c__1, (char *)&d__, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, " +/- ", (ftnlen)5);
	do_fio(&c__1, (char *)&sigd, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    return 0;
} /* doors_ */

integer fvtdo_(doublereal *chi2, integer *nt, integer *tlist)
{
    /* System generated locals */
    integer ret_val;

    /* Local variables */
    static doublereal x0[3], cx0[9]	/* was [3][3] */;
    extern integer fvfit_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    extern /* Subroutine */ int fvcopy_(doublereal *, doublereal *, integer *,
	     integer *);

/* ----------------------------------------------------------------------- */
/* -- do vertex fit. */
/* -- use tx(i,1) and tCx(i,j,1) as initial values */
/* -- and track parameters from th, tCh */
/* -- put results into tx, tCx and tq, tCq */
/* -- */
/* CC      parameter (dataDir = 'disk$user:[bauerdick.fv.dat]') */
/* -- ???? is not necessary */
    /* Parameter adjustments */
    --tlist;

    /* Function Body */
    fvcopy_(x0, fvtcom_1.tx, &c__3, &c__1);
    fvcopy_(cx0, fvtcom_1.tcx, &c__3, &c__3);
    ret_val = fvfit_(fvtcom_1.tx, fvtcom_1.tcx, fvtcom_1.tq, fvtcom_1.tcq, 
	    fvtcom_1.tchi2, chi2, nt, &tlist[1], x0, cx0, fvtcom_1.th, 
	    fvtcom_1.tgh);
    return ret_val;
} /* fvtdo_ */

integer fvtread_(integer *krun, integer *kevt)
{
    /* System generated locals */
    integer ret_val, i__1;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    , f_open(olist *), s_rsle(cilist *), do_lio(integer *, integer *, 
	    char *, ftnlen), e_rsle(void), f_clos(cllist *);

    /* Local variables */
    static integer j, i0, i1;
    static char fnam[80];
    static integer status;
    extern integer fvluinv_(doublereal *, doublereal *, integer *);

    /* Fortran I/O blocks */
    static icilist io___99 = { 0, fnam, 0, "(a,a,i5.5,a,i6.6,a)", 80, 1 };
    static cilist io___100 = { 0, 42, 0, 0, 0 };
    static cilist io___102 = { 0, 42, 0, 0, 0 };
    static cilist io___104 = { 0, 42, 0, 0, 0 };
    static cilist io___105 = { 0, 42, 0, 0, 0 };
    static cilist io___107 = { 0, 42, 0, 0, 0 };
    static cilist io___108 = { 0, 42, 0, 0, 0 };


/* ----------------------------------------------------------------------- */
/* -- read in track data, return */
/* --   tx, Cx: vertex position in carthesian coord. and cov. matrix */
/* --   th, Ch: track parameter in ALEPH comvention and cov. matrix */
/* --   tnt: number of tracks */
/* --   tGh: inverse of covariance matrix of track parameter */
/* -- in common /FVCOM/ */
/* -- */
/* -- returns ERROR if problem during inverting tCh occurred */

/* CC      parameter (dataDir = 'disk$user:[bauerdick.fv.dat]') */
    ret_val = 1;
    s_wsfi(&io___99);
    do_fio(&c__1, "dat/", (ftnlen)4);
    do_fio(&c__1, "tr", (ftnlen)2);
    do_fio(&c__1, (char *)&(*krun), (ftnlen)sizeof(integer));
    do_fio(&c__1, "e", (ftnlen)1);
    do_fio(&c__1, (char *)&(*kevt), (ftnlen)sizeof(integer));
    do_fio(&c__1, ".dat", (ftnlen)4);
    e_wsfi();
    o__1.oerr = 0;
    o__1.ounit = 42;
    o__1.ofnmlen = 80;
    o__1.ofnm = fnam;
    o__1.orl = 0;
    o__1.osta = "OLD";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    s_rsle(&io___100);
    for (i0 = 1; i0 <= 3; ++i0) {
	do_lio(&c__5, &c__1, (char *)&fvtcom_1.tx[i0 - 1], (ftnlen)sizeof(
		doublereal));
    }
    e_rsle();
    s_rsle(&io___102);
    for (i1 = 1; i1 <= 3; ++i1) {
	for (i0 = 1; i0 <= 3; ++i0) {
	    do_lio(&c__5, &c__1, (char *)&fvtcom_1.tcx[i0 + (i1 + 3) * 3 - 13]
		    , (ftnlen)sizeof(doublereal));
	}
    }
    e_rsle();
    s_rsle(&io___104);
    do_lio(&c__5, &c__1, (char *)&fvtcom_1.tw2pt, (ftnlen)sizeof(doublereal));
    e_rsle();
    s_rsle(&io___105);
    do_lio(&c__3, &c__1, (char *)&fvtcom_1.tnt, (ftnlen)sizeof(integer));
    e_rsle();
    i__1 = fvtcom_1.tnt;
    for (j = 1; j <= i__1; ++j) {
	s_rsle(&io___107);
	for (i0 = 1; i0 <= 5; ++i0) {
	    do_lio(&c__5, &c__1, (char *)&fvtcom_1.th[i0 + j * 5 - 6], (
		    ftnlen)sizeof(doublereal));
	}
	e_rsle();
	s_rsle(&io___108);
	for (i1 = 1; i1 <= 5; ++i1) {
	    for (i0 = 1; i0 <= 5; ++i0) {
		do_lio(&c__5, &c__1, (char *)&fvtcom_1.tch[i0 + (i1 + j * 5) *
			 5 - 31], (ftnlen)sizeof(doublereal));
	    }
	}
	e_rsle();
    }
    cl__1.cerr = 0;
    cl__1.cunit = 42;
    cl__1.csta = 0;
    f_clos(&cl__1);
/* -- calculate inverse of covariance matrices for th */
    i__1 = fvtcom_1.tnt;
    for (j = 1; j <= i__1; ++j) {
/* CCCCCCCCC        call fvCopy(tGh(1,1,j),tCh(1,1,j),5,5) */
	status = fvluinv_(&fvtcom_1.tgh[(j * 5 + 1) * 5 - 30], &fvtcom_1.tch[(
		j * 5 + 1) * 5 - 30], &c__5);
/* CCC        status = fvCalcG(tGh(1,1,j),tCh(1,1,j),5) */
/* CCC        if (status .EQ. INFO) then */
/* CCC          write(fvtPlun,'(1x,a)') */
/* CCC     1      'editing of singular values during inv of Ch0' */
/* CCC          call fvAB(one,tGh(1,1,j),tCh(1,1,j),5,5,5) */
/* CCC          write(fvtPlun,'(/1x,5g10.3)') ((one(i0,i1), i0=1,5),i1=1,5) */
/* CCC        end if */
	if ((status & 1) != 1) {
	    ret_val = status;
	}
    }
    return ret_val;
} /* fvtread_ */

integer fvtinvm_(doublereal *mi, doublereal *sigmi, integer *nt, integer *
	tlist, doublereal *t, doublereal *ct)
{
    /* System generated locals */
    integer ret_val, i__1;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, i0, i1;
    static doublereal pi[4];
    static integer it;
    static doublereal cpi[16]	/* was [4][4] */;
    extern /* Subroutine */ int fvzeroa_(doublereal *, integer *, integer *);

/* ----------------------------------------------------------------------- */
/* -- calculate invariant mass and error */
/* -- from 4-momenta t, Ct */

/* CC      parameter (dataDir = 'disk$user:[bauerdick.fv.dat]') */
    /* Parameter adjustments */
    --tlist;
    t -= 5;
    ct -= 21;

    /* Function Body */
    fvzeroa_(pi, &c__4, &c__1);
    fvzeroa_(cpi, &c__4, &c__4);
    i__1 = *nt;
    for (i__ = 1; i__ <= i__1; ++i__) {
	it = tlist[i__];
	for (i0 = 1; i0 <= 4; ++i0) {
	    pi[i0 - 1] += t[i0 + (it << 2)];
	    for (i1 = 1; i1 <= 4; ++i1) {
		cpi[i0 + (i1 << 2) - 5] += ct[i0 + (i1 + (it << 2) << 2)];
	    }
	}
    }
/* Computing 2nd power */
    d__1 = pi[3];
/* Computing 2nd power */
    d__2 = pi[0];
/* Computing 2nd power */
    d__3 = pi[1];
/* Computing 2nd power */
    d__4 = pi[2];
    *mi = sqrt(d__1 * d__1 - d__2 * d__2 - d__3 * d__3 - d__4 * d__4);
    *sigmi = pi[0] * cpi[0] * pi[0] + pi[1] * cpi[5] * pi[1] + pi[2] * cpi[10]
	     * pi[2] + pi[3] * cpi[15] * pi[3] + (pi[0] * (cpi[4] * pi[1] + 
	    cpi[8] * pi[2] - cpi[12] * pi[3]) + pi[1] * (cpi[9] * pi[2] - cpi[
	    13] * pi[3]) - pi[2] * cpi[14] * pi[3]) * 2.f;
    *sigmi = sqrt(*sigmi) / *mi;
    return ret_val;
} /* fvtinvm_ */

/* Subroutine */ int fvthelix2p4_(doublereal *t, doublereal *ct, integer *nt, 
	integer *tlist, integer *plist, doublereal *p, doublereal *ch)
{
    /* Initialized data */

    static doublereal mass[4] = { 5.1099906e-4,.105658387,1.7841,.1395675 };

    /* System generated locals */
    integer i__1;

    /* Local variables */
    extern /* Subroutine */ int fvenergy_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    static integer i__;
    static doublereal m;
    extern /* Subroutine */ int fvhelix2p4_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    static integer ip, it;

/* ----------------------------------------------------------------------- */
/* -- calculate 4-momentum from track parameters */
/* -- using track/particle list tList, pList */

/* CC      parameter (dataDir = 'disk$user:[bauerdick.fv.dat]') */
    /* Parameter adjustments */
    t -= 5;
    ct -= 21;
    --plist;
    --tlist;
    p -= 6;
    ch -= 31;

    /* Function Body */
    i__1 = *nt;
    for (i__ = 1; i__ <= i__1; ++i__) {
	it = tlist[i__];
	fvhelix2p4_(&t[(it << 2) + 1], &ct[((it << 2) + 1 << 2) + 1], &
		fvtcom_1.tw2pt, &p[it * 5 + 1], &ch[(it * 5 + 1) * 5 + 1]);
	ip = plist[i__];
	m = mass[ip - 1];
	fvenergy_(&t[(it << 2) + 1], &ct[((it << 2) + 1 << 2) + 1], &m, &t[(
		it << 2) + 1], &ct[((it << 2) + 1 << 2) + 1]);
    }
    return 0;
} /* fvthelix2p4_ */

/* Subroutine */ int fvtqm2p4_(doublereal *tl, doublereal *ctl, integer *nt, 
	integer *tlist, integer *plist, doublereal *ql, doublereal *cql)
{
    /* Initialized data */

    static doublereal mass[4] = { 5.1099906e-4,.105658387,1.7841,.1395675 };

    /* System generated locals */
    integer i__1;

    /* Local variables */
    extern /* Subroutine */ int fvenergy_();
    static integer i__;
    static doublereal m;
    static integer ip, it;
    extern /* Subroutine */ int fvq2p4_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);

/* ----------------------------------------------------------------------- */
/* -- calculate list of 4-momentum vectors */
/* -- from mass and q-vectors {w, tl, psi} */
/* -- using track/particle list tList, pList */

/* CC      parameter (dataDir = 'disk$user:[bauerdick.fv.dat]') */
    /* Parameter adjustments */
    tl -= 5;
    ctl -= 21;
    --plist;
    --tlist;
    ql -= 4;
    cql -= 13;

    /* Function Body */
    i__1 = *nt;
    for (i__ = 1; i__ <= i__1; ++i__) {
	it = tlist[i__];
	fvq2p4_(&tl[(it << 2) + 1], &ctl[((it << 2) + 1 << 2) + 1], &ql[it * 
		3 + 1], &cql[(it * 3 + 1) * 3 + 1], &fvtcom_1.tw2pt);
	ip = plist[i__];
	m = mass[ip - 1];
	fvenergy_(&tl[(it << 2) + 1], &ctl[((it << 2) + 1 << 2) + 1], &m);
    }
    return 0;
} /* fvtqm2p4_ */

integer fvtvert_(doublereal *v, doublereal *cv, doublereal *q, doublereal *cq,
	 doublereal *v0, doublereal *cv0, doublereal *q0, doublereal *cq0)
{
    /* System generated locals */
    integer ret_val;

/* ----------------------------------------------------------------------- */
/* -- */
/* CC      parameter (dataDir = 'disk$user:[bauerdick.fv.dat]') */
    /* Parameter adjustments */
    cq0 -= 4;
    --q0;
    cv0 -= 4;
    --v0;
    cq -= 4;
    --q;
    cv -= 4;
    --v;

    /* Function Body */
    ret_val = 1;
    return ret_val;
} /* fvtvert_ */

integer fvttav_(doublereal *x, doublereal *cx, doublereal *p, doublereal *cp, 
	doublereal *chi2, doublereal *x0, doublereal *cx0, doublereal *h__, 
	doublereal *ch)
{
    /* System generated locals */
    integer ret_val;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    double sqrt(doublereal);

    /* Local variables */
    extern integer fvfilter_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    static doublereal a[15]	/* was [5][3] */, b[15]	/* was [5][3] */, c__[
	    9]	/* was [3][3] */, d__[9]	/* was [3][3] */, e[9]	/* 
	    was [3][3] */;
    static integer i__;
    static doublereal q[3], v[3], c0[9]	/* was [3][3] */, d0[9]	/* was [3][3] 
	    */, h0[5], q0[3], v0[3], gv[9]	/* was [3][3] */;
    extern integer fvxp_(doublereal *, doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *), fvv0q0_(doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *), fvabh0_(doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    static integer status;

    /* Fortran I/O blocks */
    static cilist io___140 = { 0, 0, 0, "(1x,g10.3)", 0 };
    static cilist io___141 = { 0, 0, 0, "(1x,3(g10.3,a,g10.3))", 0 };
    static cilist io___143 = { 0, 0, 0, "(1x,3(g10.3,a,g10.3))", 0 };
    static cilist io___144 = { 0, 0, 0, "(1x,3(g10.3,a,g10.3))", 0 };
    static cilist io___145 = { 0, 0, 0, "(1x,3(g10.3,a,g10.3))", 0 };
    static cilist io___146 = { 0, 0, 0, "(1x)", 0 };


/* ----------------------------------------------------------------------- */
/* -- track and vertex: */
/* -- do the refit of track with */
/* -- helix parameters h and covariance matrix Ch, */
/* -- using point x0 and its covariance matrix Cx0 */
/* -- h = {w, tl, psi0, d0, z0} and x0 = {x, y, z} */
/* -- return results in x, Cx, and p, Cp */
/* -- */
/* CC      parameter (dataDir = 'disk$user:[bauerdick.fv.dat]') */
    /* Parameter adjustments */
    ch -= 6;
    --h__;
    cx0 -= 4;
    --x0;
    cp -= 6;
    --p;
    cx -= 4;
    --x;

    /* Function Body */
    ret_val = 1;
/* -- calculate vertex and 3-momentum */
/* -- at approximate vertex position x0 */
/* -- and the appropriate covariance matrices */
    status = fvv0q0_(v0, c0, q0, d0, &x0[1], &cx0[4], &h__[1], &ch[6]);
    if ((status & 1) != 1) {
	ret_val = status;
	return ret_val;
    }
/* -- calculate coefficients for measurement equation, A, B, h0 */
    status = fvabh0_(a, b, h0, v0, q0);
    if ((status & 1) != 1) {
	ret_val = status;
	return ret_val;
    }
/* -- calculate filtered state vectors and covariance matrices */
    status = fvfilter_(v, c__, gv, q, d__, e, chi2, v0, c0, &h__[1], &ch[6], 
	    a, b, h0);
    if ((status & 1) != 1) {
	ret_val = status;
	return ret_val;
    }
/* -- re-calculate coefficients for measurement equation, A, B, h0 */
/* -- and check on difference, ev. re-iterate (not done at the moment) */
    status = fvabh0_(a, b, h0, v, q);
    if ((status & 1) != 1) {
	ret_val = status;
	return ret_val;
    }
/* -- calculate estimated vertex and track parameters */
    status = fvxp_(&x[1], &cx[4], &p[1], &cp[6], v, q, c__, d__, e, a, b);
/* -- some printout */
    io___140.ciunit = fvtcom_1.fvtplun;
    s_wsfe(&io___140);
    do_fio(&c__1, (char *)&(*chi2), (ftnlen)sizeof(doublereal));
    e_wsfe();
    io___141.ciunit = fvtcom_1.fvtplun;
    s_wsfe(&io___141);
    for (i__ = 1; i__ <= 3; ++i__) {
	do_fio(&c__1, (char *)&v0[i__ - 1], (ftnlen)sizeof(doublereal));
	do_fio(&c__1, " +/-", (ftnlen)4);
	d__1 = sqrt(c0[i__ + i__ * 3 - 4]);
	do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
    }
    e_wsfe();
    io___143.ciunit = fvtcom_1.fvtplun;
    s_wsfe(&io___143);
    for (i__ = 1; i__ <= 3; ++i__) {
	do_fio(&c__1, (char *)&v[i__ - 1], (ftnlen)sizeof(doublereal));
	do_fio(&c__1, " +/-", (ftnlen)4);
	d__1 = sqrt(c__[i__ + i__ * 3 - 4]);
	do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
    }
    e_wsfe();
/* CC      write(fvtPlun,'(1x,3(g10.3,a,g10.3))') (x0(i), ' +/-', */
/* CC     1   sqrt(Cx0(i,i)) ,i=1,3) */
/* CC      write(fvtPlun,'(1x,3(g10.3,a,g10.3))') (x(i), ' +/-', */
/* CC     1   sqrt(Cx(i,i)) ,i=1,3) */
    io___144.ciunit = fvtcom_1.fvtplun;
    s_wsfe(&io___144);
    for (i__ = 1; i__ <= 3; ++i__) {
	do_fio(&c__1, (char *)&q0[i__ - 1], (ftnlen)sizeof(doublereal));
	do_fio(&c__1, " +/-", (ftnlen)4);
	d__1 = sqrt(ch[i__ + i__ * 5]);
	do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
    }
    e_wsfe();
    io___145.ciunit = fvtcom_1.fvtplun;
    s_wsfe(&io___145);
    for (i__ = 1; i__ <= 3; ++i__) {
	do_fio(&c__1, (char *)&q[i__ - 1], (ftnlen)sizeof(doublereal));
	do_fio(&c__1, " +/-", (ftnlen)4);
	d__1 = sqrt(d__[i__ + i__ * 3 - 4]);
	do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
    }
    e_wsfe();
    io___146.ciunit = fvtcom_1.fvtplun;
    s_wsfe(&io___146);
    e_wsfe();
    return ret_val;
} /* fvttav_ */

/* Subroutine */ int errinit_(void)
{
    return 0;
} /* errinit_ */

/* Subroutine */ int errsumm_(void)
{
    return 0;
} /* errsumm_ */

/* Main program alias */ int fvt_ () { MAIN__ (); return 0; }
