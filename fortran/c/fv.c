/* fv.f -- translated by f2c (version 20090411).
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

struct fvcom_1_ {
    logical print;
    integer plun;
    char text[132];
};

#define fvcom_1 (*(struct fvcom_1_ *) &fvcom_)

/* Initialized data */

struct {
    logical e_1;
    integer e_2;
    integer fill_3[33];
    } fvcom_ = { FALSE_, 41 };


/* Table of constant values */

static integer c__3 = 3;
static integer c__1 = 1;
static integer c__2 = 2;
static integer c__5 = 5;
static doublereal c_b86 = 0.;
static integer c__4 = 4;
static doublereal c_b134 = 1.;
static doublereal c_b509 = 6.28318530717959;
static integer c__6 = 6;
static integer c__20 = 20;
static integer c__9 = 9;

/* ---------------------------------------------------------------------- */

/* FV, package to fit vertices and find secondary vertices */

/* author: Lothar A.T. Bauerdick, CERN/PPE */
/* date:   Jan. 1991 */

/* ref.:   R. Fruehwirt, ``Applications of Filter Methods...'' */
/*         HEPHY-PUB 516/88, Vienna */


/* some notation: */

/* w:  omega, 1/R curvature >0 if helix turns anti-clockwise */
/* tl: tan(lambda), tangens of dip angle */
/* psi: azimuth angle (in r-phi plane) */
/* psi0: psi @ point of closest approach to coordinate origin */
/* d0:   distance of helix at point of closest approach to origin in r-phi plane */
/*       sign convention of angular momentum Lz */
/* z0:   z coordinate of point of closest approach of helix to origin */

/* ---------------------------------------------------------------------- */
/* Subroutine */ int fv_(void)
{
    return 0;
} /* fv_ */

/* Subroutine */ int fvblkd_(void)
{
    return 0;
} /* fvblkd_ */


integer fvfit_(doublereal *x, doublereal *cx, doublereal *ql, doublereal *cql,
	 doublereal *chi2l, doublereal *chi2t, integer *nt, integer *tlist, 
	doublereal *x0, doublereal *cx0, doublereal *hl, doublereal *ghl)
{
    /* System generated locals */
    integer ret_val, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    double sqrt(doublereal);
    integer s_wsfi(icilist *), e_wsfi(void);

    /* Local variables */
    extern integer fvfilter_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), 
	    fvremove_(doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), 
	    fvsmooth_(doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    static doublereal c__[9]	/* was [3][3] */, e[9]	/* was [3][3] */;
    static integer i__;
    static doublereal v[3], c0[9]	/* was [3][3] */;
    static integer i0;
    static doublereal v0[3];
    static integer it;
    static doublereal gv[9]	/* was [3][3] */, pp[5], gv0[9]	/* was [3][3] 
	    */, cpp[25]	/* was [5][5] */, chi2;
    extern integer fvhch_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    extern doublereal fvprob_(doublereal *, integer *);
    extern /* Subroutine */ int fvcopy_(doublereal *, doublereal *, integer *,
	     integer *);
    static integer status;
    extern integer fvcalcg_(doublereal *, doublereal *, integer *);
    extern /* Subroutine */ int fverror_(char *, ftnlen);
    extern integer fvluinv_(doublereal *, doublereal *, integer *);

    /* Fortran I/O blocks */
    static cilist io___3 = { 0, 0, 0, "(1x,a)", 0 };
    static cilist io___4 = { 0, 0, 0, "(1x,3(g10.3,a,g10.3))", 0 };
    static icilist io___15 = { 0, fvcom_1.text, 0, "(a,i4)", 132, 1 };
    static cilist io___16 = { 0, 0, 0, "(1x,a,i5)", 0 };
    static icilist io___17 = { 0, fvcom_1.text, 0, "(a,i4)", 132, 1 };
    static cilist io___18 = { 0, 0, 0, "(1x,a,g10.3,a,g10.3)", 0 };
    static cilist io___20 = { 0, 0, 0, "(1x,5(g10.3,a,g10.3))", 0 };
    static cilist io___22 = { 0, 0, 0, "(1x,5(g10.3,a,g10.3))", 0 };
    static cilist io___23 = { 0, 0, 0, "(1x,a)", 0 };
    static cilist io___24 = { 0, 0, 0, "(1x,a,g10.3,a,i3,a,i3,a,g10.3,a,g10."
	    "3,a)", 0 };
    static cilist io___25 = { 0, 0, 0, "(1x,3(g10.3,a,g10.3))", 0 };


/* ----------------------------------------------------------------------- */
/* -- calculate vertex x and its covariance matrix Cx (x = {x,y,z}), */
/* -- the list of momentum vectors ql and their covariance matrices */
/* -- Cql (q = {w, tl, psi}) and the total chi-square chi2t */
/* -- for nt tracks with numbers in tList */
/* -- from initial value for vertex x0, Cx0 and */
/* -- from list of track parameters hl and list of */
/* -- inverse of covariance matrices of track parameters Ghl */
/* -- 3-momenta q are returned in {w, tl, psi}-system */
/* --   (we do not know the transformation of w->pt at this point!) */
/* -- */
/* -- input: */
/* --     nt             : number of tracks in tList */
/* --     tList({1..nt}) : index of tracks in hl, Ghl */
/* --     x0({1..3})     : initial guess of vertex position x0{x,y,z} */
/* --     Cx0({1..3},{1..3}) : covariance matrix of x0 */
/* --     hl({1..5}, tList(i)) : track parameters */
/* --                      p0={w,tl,psi0,d0,z0} of track i */
/* --     Ghl({1..5},{1..5}, tList(i)) : inverse of covariance matrix */
/* --                      (cov(p0(i)))^(-1) of p0 of track i */
/* -- output: */
/* --     x({1..3})              <- v({x,y,z}), coords of fitted vertx */
/* --     Cx({1..3})             <- cov{v} covariance matrix of v */
/* --     ql({1..3}, tList(i))   <- q({w,tl,psi}) */
/* --     Cql({1..3},{1..3},tList(i)) <- covariance matrix of q */
/* --     chi2l(tList(i))        <- chi2 for this track belonging to v */
/* --     chi2t                  <- total chi2 for the vertex fit */
/* -- */
/* -- author: Lothar A.T. Bauerdick, CERN/PPE */
/* -- date:   Jan. 1991 */
/* -- */
    /* Parameter adjustments */
    --x;
    cx -= 4;
    ql -= 4;
    cql -= 13;
    --chi2l;
    --tlist;
    --x0;
    cx0 -= 4;
    hl -= 6;
    ghl -= 31;

    /* Function Body */
    ret_val = 1;

/* -- the filter */

/* -- calculate start values v0, C0 from x0, Cx0 */
/* -- v, Cv are in carth. coords too, ergo just copy */
    fvcopy_(v0, &x0[1], &c__3, &c__1);
    fvcopy_(c0, &cx0[4], &c__3, &c__3);
    if (fvcom_1.print) {
	io___3.ciunit = fvcom_1.plun;
	s_wsfe(&io___3);
	do_fio(&c__1, "start values for Filter", (ftnlen)23);
	e_wsfe();
	io___4.ciunit = fvcom_1.plun;
	s_wsfe(&io___4);
	for (i0 = 1; i0 <= 3; ++i0) {
	    do_fio(&c__1, (char *)&v0[i0 - 1], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, " +/-", (ftnlen)4);
	    d__1 = sqrt(c0[i0 + i0 * 3 - 4]);
	    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	}
	e_wsfe();
    }
/* -- calculate inverse of covariance matrices for v0 */
    status = fvcalcg_(gv0, c0, &c__3);
    if ((status & 1) != 1) {
	ret_val = status;
	return ret_val;
    }
    i__1 = *nt;
    for (i__ = 1; i__ <= i__1; ++i__) {
	it = tlist[i__];
	if (it != 0) {
	    status = fvfilter_(v, c__, gv, &ql[it * 3 + 1], &cql[(it * 3 + 1) 
		    * 3 + 1], e, &chi2, v0, gv0, &hl[it * 5 + 1], &ghl[(it * 
		    5 + 1) * 5 + 1]);
	    if ((status & 1) != 1) {
		s_wsfi(&io___15);
		do_fio(&c__1, "Problems in fvFilter for Track ", (ftnlen)31);
		do_fio(&c__1, (char *)&it, (ftnlen)sizeof(integer));
		e_wsfi();
		fverror_(fvcom_1.text, (ftnlen)132);
		ret_val = status;
	    }
/* -- use v and Gv as input vertex v0, Gv0 for next track */
	    fvcopy_(v0, v, &c__3, &c__1);
	    fvcopy_(gv0, gv, &c__3, &c__3);
	}
    }

/* -- the smoother */
/* -- v, C, Gv contain vertex, cov. matrix and cov^(-1) */
    *chi2t = 0.;
    i__1 = *nt;
    for (i__ = 1; i__ <= i__1; ++i__) {
	it = tlist[i__];
	if (it != 0) {
/* -- some printout */
	    if (fvcom_1.print) {
		io___16.ciunit = fvcom_1.plun;
		s_wsfe(&io___16);
		do_fio(&c__1, "Smoother for Track ", (ftnlen)19);
		do_fio(&c__1, (char *)&it, (ftnlen)sizeof(integer));
		e_wsfe();
	    }
	    status = fvsmooth_(&ql[it * 3 + 1], &cql[(it * 3 + 1) * 3 + 1], e,
		     &chi2, v, c__, &ql[it * 3 + 1], &hl[it * 5 + 1], &ghl[(
		    it * 5 + 1) * 5 + 1]);
	    if ((status & 1) != 1) {
		s_wsfi(&io___17);
		do_fio(&c__1, "Problems in smoother for track ", (ftnlen)31);
		do_fio(&c__1, (char *)&it, (ftnlen)sizeof(integer));
		e_wsfi();
		fverror_(fvcom_1.text, (ftnlen)132);
		ret_val = status;
	    }
/* -- calculate total chi2 */
	    *chi2t += chi2;
	    if (fvcom_1.print) {
/* -- print out some results */
		io___18.ciunit = fvcom_1.plun;
		s_wsfe(&io___18);
		do_fio(&c__1, "chi2 =", (ftnlen)6);
		do_fio(&c__1, (char *)&chi2, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, "prob (2 d.o.f.) =", (ftnlen)17);
		d__1 = fvprob_(&chi2, &c__2);
		do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
		e_wsfe();
/* -- get old track parameters error matrix */
/* CCCCCCCCC            call fvCopy(Cpp,Ghl(1,1,it),dh,dh) */
		status = fvluinv_(cpp, &ghl[(it * 5 + 1) * 5 + 1], &c__5);
		if ((status & 1) != 1) {
		    ret_val = status;
		}
		io___20.ciunit = fvcom_1.plun;
		s_wsfe(&io___20);
		for (i0 = 1; i0 <= 5; ++i0) {
		    do_fio(&c__1, (char *)&hl[i0 + it * 5], (ftnlen)sizeof(
			    doublereal));
		    do_fio(&c__1, " +/-", (ftnlen)4);
		    d__1 = sqrt(cpp[i0 + i0 * 5 - 6]);
		    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
		}
		e_wsfe();
/* -- calculate new track parameters and error matrix for printout */
		status = fvhch_(pp, cpp, v, &ql[it * 3 + 1], c__, &cql[(it * 
			3 + 1) * 3 + 1], e);
		if ((status & 1) != 1) {
		    ret_val = status;
		}
		io___22.ciunit = fvcom_1.plun;
		s_wsfe(&io___22);
		for (i0 = 1; i0 <= 5; ++i0) {
		    do_fio(&c__1, (char *)&pp[i0 - 1], (ftnlen)sizeof(
			    doublereal));
		    do_fio(&c__1, " +/-", (ftnlen)4);
		    d__1 = sqrt(cpp[i0 + i0 * 5 - 6]);
		    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
		}
		e_wsfe();
	    }
	}
    }
/* -- save smoothed vertex position in x = {x, y, z} */
    fvcopy_(&x[1], v, &c__3, &c__1);
    fvcopy_(&cx[4], c__, &c__3, &c__3);
/* -- smoothed vertex position still in v, C, Gv */
    if (fvcom_1.print) {
	io___23.ciunit = fvcom_1.plun;
	s_wsfe(&io___23);
	do_fio(&c__1, "Smoothed result", (ftnlen)15);
	e_wsfe();
	io___24.ciunit = fvcom_1.plun;
	s_wsfe(&io___24);
	do_fio(&c__1, "chi2 =", (ftnlen)6);
	do_fio(&c__1, (char *)&(*chi2t), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, "prob for ", (ftnlen)9);
	i__1 = *nt << 1;
	do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
	do_fio(&c__1, " (", (ftnlen)2);
	i__2 = (*nt << 1) - 3;
	do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(integer));
	do_fio(&c__1, ") d.o.f. =", (ftnlen)10);
	i__3 = *nt << 1;
	d__1 = fvprob_(chi2t, &i__3);
	do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, " (", (ftnlen)2);
	i__4 = (*nt << 1) - 3;
	d__2 = fvprob_(chi2t, &i__4);
	do_fio(&c__1, (char *)&d__2, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, ")", (ftnlen)1);
	e_wsfe();
	io___25.ciunit = fvcom_1.plun;
	s_wsfe(&io___25);
	for (i0 = 1; i0 <= 3; ++i0) {
	    do_fio(&c__1, (char *)&v[i0 - 1], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, " +/-", (ftnlen)4);
	    d__1 = sqrt(c__[i0 + i0 * 3 - 4]);
	    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	}
	e_wsfe();
    }
/* -- check for outliers */
/* -- calculate chi2 for each track to belong to this vertex */
/* -- set alpha = cut on probaility to zero, do not remove any track */
    if (*nt <= 2) {
	chi2l[tlist[1]] = 0.;
	chi2l[tlist[2]] = 0.;
    } else {
	status = fvremove_(v0, c0, gv0, &ql[4], &cql[13], &chi2l[1], nt, &
		tlist[1], &c_b86, v, gv, &hl[6], &ghl[31]);
	if ((status & 1) != 1) {
	    ret_val = status;
	}
    }
    return ret_val;
} /* fvfit_ */

integer fvfitd_(doublereal *tau, doublereal *sigtau, doublereal *chi2, 
	doublereal *tau0, doublereal *v0, doublereal *c0, doublereal *vp, 
	doublereal *cp, doublereal *p, doublereal *dp, doublereal *ep)
{
    /* System generated locals */
    integer ret_val;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    extern integer fvsvdinv_(doublereal *, integer *);
    static doublereal a[3], c__[9]	/* was [3][3] */, g[9]	/* was [3][3] 
	    */, v[3], jp[9]	/* was [3][3] */, pp, jvp[9]	/* was [3][3] 
	    */;
    static integer iter;
    static doublereal taup;
    extern /* Subroutine */ int fvapb_(doublereal *, doublereal *, doublereal 
	    *, integer *, integer *);
    static doublereal temp11[1], temp33[9]	/* was [3][3] */;
    extern /* Subroutine */ int fvabct_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *), fvatbc_(
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *), fvapbt_(doublereal *, doublereal *, doublereal *, 
	    integer *, integer *);
    static integer status;
    static doublereal chi2old;
    extern /* Subroutine */ int fvsabat_(doublereal *, doublereal *, 
	    doublereal *, integer *, integer *), fvsatba_(doublereal *, 
	    doublereal *, doublereal *, integer *, integer *), fvjfitd_(
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);

    /* Fortran I/O blocks */
    static cilist io___36 = { 0, 0, 0, "(1x,a)", 0 };


/* ----------------------------------------------------------------------- */
/* -- */
/* -- fit tau for */
/* --    v = vp - tau*p/|p| */
/* -- where v is measured as v0 with covariance C0 */
/* -- and vp,p are measured with Cp=cov{vp,vp},Dp=cov{p,p},Ep=cov{vp,p} */
    /* Parameter adjustments */
    ep -= 4;
    dp -= 4;
    --p;
    cp -= 4;
    --vp;
    c0 -= 4;
    --v0;

    /* Function Body */
    ret_val = 1;
    *sigtau = 0.;
/* Computing 2nd power */
    d__1 = p[1];
/* Computing 2nd power */
    d__2 = p[2];
/* Computing 2nd power */
    d__3 = p[3];
    pp = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
/* -- A = dv / dtau */
    a[0] = -p[1] / pp;
    a[1] = -p[2] / pp;
    a[2] = -p[3] / pp;
    *tau = *tau0;
    iter = 0;
    chi2old = 1e10;
L1:
/* -- calculate jacobian J(v|vp,p) and transform covariance matrix */
/* -- for measurement of secondary vertex vp/p for given tau */
    fvjfitd_(jvp, jp, &p[1], &pp, tau);
    fvsabat_(c__, jvp, &cp[4], &c__3, &c__3);
    fvsabat_(temp33, jp, &dp[4], &c__3, &c__3);
    fvapb_(c__, c__, temp33, &c__3, &c__3);
    fvabct_(temp33, jvp, &ep[4], jp, &c__3, &c__3);
    fvapb_(c__, c__, temp33, &c__3, &c__3);
    fvapbt_(c__, c__, temp33, &c__3, &c__3);
/* -- total covariance matrix is C + C0, invert it */
    fvapb_(g, &c0[4], c__, &c__3, &c__3);
    status = fvsvdinv_(g, &c__3);
    if (status == 3) {
	io___36.ciunit = fvcom_1.plun;
	s_wsfe(&io___36);
	do_fio(&c__1, "editing of singular values done during inv. of in tau"
		" fit", (ftnlen)57);
	e_wsfe();
    }
    if ((status & 1) != 1) {
	ret_val = status;
	return ret_val;
    }
/* -- calculate chi2 */
    v[0] = v0[1] - (vp[1] - *tau * p[1] / pp);
    v[1] = v0[2] - (vp[2] - *tau * p[2] / pp);
    v[2] = v0[3] - (vp[3] - *tau * p[3] / pp);
    fvsatba_(temp11, v, g, &c__3, &c__1);
    *chi2 = temp11[0];
    ++iter;
    if ((d__1 = chi2old - *chi2, abs(d__1)) < .1) {
/* -- chi^2 does not change anymore, get out of here */
	goto L99;
    }
    if (iter > 100) {
/* -- doesn't converge, return ERROR */
	ret_val = 2;
	return ret_val;
    }
    chi2old = *chi2;
/* -- do the fit */
/* -- tau(i) = tau(i-1) + {A^T.G.A}^(-1).A^T.G.(v0-v(tau(i-1))) */
/* -- where A = dv/dtau | tau=tau(i-1) */
/* -- and G = (cov{v})^(-1) */
    fvsatba_(temp11, a, g, &c__3, &c__1);
    taup = 1. / temp11[0];
    fvatbc_(temp11, a, g, v, &c__3, &c__1);
    taup = *tau + taup * temp11[0];
    *tau = taup;
    goto L1;
L99:
/* -- re-calculate cov with most recent G */
    fvsatba_(temp11, a, g, &c__3, &c__1);
/* Computing MAX */
    d__1 = 0., d__2 = 1. / temp11[0];
    *sigtau = sqrt((max(d__1,d__2)));
    return ret_val;
} /* fvfitd_ */

/* Subroutine */ int fvjfitd_(doublereal *jvp, doublereal *jp, doublereal *p, 
	doublereal *pp, doublereal *tau)
{
    /* System generated locals */
    doublereal d__1;

    /* Local variables */
    static doublereal pp2, tpp;

/* ----------------------------------------------------------------------- */
/* -- calculate Jacobian for v = vp - tau*p/|p| */
/* -- Jvp = J(v|vp), Jp = J(v|p) */
/* -- */
    /* Parameter adjustments */
    --p;
    jp -= 4;
    jvp -= 4;

    /* Function Body */
/* Computing 2nd power */
    d__1 = *pp;
    pp2 = d__1 * d__1;
    tpp = *tau / *pp;
    jvp[4] = 1.;
    jvp[7] = 0.;
    jvp[10] = 0.;
    jvp[5] = 0.;
    jvp[8] = 1.;
    jvp[11] = 0.;
    jvp[6] = 0.;
    jvp[9] = 0.;
    jvp[12] = 1.;
    jp[4] = -tpp * (1. - p[1] * p[1] / pp2);
    jp[7] = tpp * p[1] * p[2] / pp2;
    jp[10] = tpp * p[1] * p[3] / pp2;
    jp[5] = jp[7];
    jp[8] = -tpp * (1. - p[2] * p[2] / pp2);
    jp[11] = tpp * p[2] * p[3] / pp2;
    jp[6] = jp[10];
    jp[9] = jp[11];
    jp[12] = -tpp * (1. - p[3] * p[3] / pp2);
    return 0;
} /* fvjfitd_ */

integer fvdist_(doublereal *d__, doublereal *sigd, doublereal *v0, doublereal 
	*cv0, doublereal *v1, doublereal *cv1)
{
    /* System generated locals */
    integer ret_val;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal dd[3], temp11[1];
    extern /* Subroutine */ int fvatbc_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *);

/* ----------------------------------------------------------------------- */
/* -- calculate distance between two vertices */

    /* Parameter adjustments */
    --cv1;
    --v1;
    cv0 -= 4;
    --v0;

    /* Function Body */
    ret_val = 1;
/* -- distance */
/* Computing 2nd power */
    d__1 = v0[1] - v1[1];
/* Computing 2nd power */
    d__2 = v0[2] - v1[2];
/* Computing 2nd power */
    d__3 = v0[3] - v1[3];
    *d__ = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
/* -- error */
    dd[0] = (v0[1] - v1[1]) / *d__;
    dd[1] = (v0[2] - v1[2]) / *d__;
    dd[2] = (v0[3] - v1[3]) / *d__;
    fvatbc_(temp11, dd, &cv0[4], dd, &c__3, &c__1);
    *sigd = temp11[0];
    fvatbc_(temp11, dd, &cv1[1], dd, &c__3, &c__1);
    *sigd = sqrt(*sigd + temp11[0]);
    return ret_val;
} /* fvdist_ */

integer fvinvm_(doublereal *m, doublereal *sigm, integer *nt, integer *tlist, 
	doublereal *w2pt, doublereal *ml, doublereal *v, doublereal *cv, 
	doublereal *ql, doublereal *cql, doublereal *ghl)
{
    /* System generated locals */
    integer ret_val, i__1;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    extern /* Subroutine */ int fvenergy_(doublereal *, doublereal *, 
	    doublereal *);
    static integer i__;
    static doublereal p[4];
    static integer i0, i1;
    static doublereal cp[16]	/* was [4][4] */, pi[4];
    static integer it;
    static doublereal cpi[16]	/* was [4][4] */;
    extern /* Subroutine */ int fvq2p4_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *), fvzeroa_(doublereal *, 
	    integer *, integer *);

/* ----------------------------------------------------------------------- */
/* -- return invariant mass of tracks */

    /* Parameter adjustments */
    --tlist;
    --ml;
    --v;
    cv -= 4;
    ql -= 4;
    cql -= 13;
    ghl -= 31;

    /* Function Body */
    ret_val = 1;
    fvzeroa_(pi, &c__4, &c__1);
    fvzeroa_(cpi, &c__4, &c__4);
    i__1 = *nt;
    for (i__ = 1; i__ <= i__1; ++i__) {
	it = tlist[i__];
	fvq2p4_(p, cp, &ql[it * 3 + 1], &cql[(it * 3 + 1) * 3 + 1], w2pt);
	fvenergy_(p, cp, &ml[it]);
	for (i0 = 1; i0 <= 4; ++i0) {
	    pi[i0 - 1] += p[i0 - 1];
	    for (i1 = 1; i1 <= 4; ++i1) {
		cpi[i0 + (i1 << 2) - 5] += cp[i0 + (i1 << 2) - 5];
	    }
	}
    }
/* -- this does not take correlations between tracks into account */
/* Computing 2nd power */
    d__1 = pi[3];
/* Computing 2nd power */
    d__2 = pi[0];
/* Computing 2nd power */
    d__3 = pi[1];
/* Computing 2nd power */
    d__4 = pi[2];
    *m = d__1 * d__1 - d__2 * d__2 - d__3 * d__3 - d__4 * d__4;
    *m = sqrt((max(*m,0.)));
    if (*m > 0.) {
	*sigm = pi[0] * cpi[0] * pi[0] + pi[1] * cpi[5] * pi[1] + pi[2] * cpi[
		10] * pi[2] + pi[3] * cpi[15] * pi[3] + (pi[0] * (cpi[4] * pi[
		1] + cpi[8] * pi[2] - cpi[12] * pi[3]) + pi[1] * (cpi[9] * pi[
		2] - cpi[13] * pi[3]) - pi[2] * cpi[14] * pi[3]) * 2.f;
	*sigm = sqrt((max(*sigm,0.))) / *m;
    } else {
	*sigm = 0.;
    }
    return ret_val;
} /* fvinvm_ */

integer fvsumq_(doublereal *p, doublereal *dp, doublereal *ep, integer *new__,
	 integer *nt, integer *tlist, doublereal *w2pt, doublereal *v, 
	doublereal *c__, doublereal *ql, doublereal *cql, doublereal *hl, 
	doublereal *ghl)
{
    /* System generated locals */
    integer ret_val, i__1, i__2;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *);

    /* Local variables */
    static doublereal d__[9]	/* was [3][3] */, e[9]	/* was [3][3] */, q[3]
	    ;
    static integer i0;
    static doublereal j0[9]	/* was [3][3] */, j1[9]	/* was [3][3] */;
    static integer i1;
    static doublereal ch[25]	/* was [5][5] */;
    extern integer fvjacobq2p_(doublereal *, doublereal *, doublereal *);
    static integer it0, it1;
    extern /* Subroutine */ int fvp2q_(doublereal *, doublereal *, doublereal 
	    *, doublereal *);
    extern integer fvhch_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    extern /* Subroutine */ int fvapb_(doublereal *, doublereal *, doublereal 
	    *, integer *, integer *), fvabt_(doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *);
    static doublereal temp31[3]	/* was [3][1] */, temp33[9]	/* was [3][3] 
	    */;
    extern /* Subroutine */ int fvq2p3_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    static doublereal charge;
    extern /* Subroutine */ int fvabct_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *), fvapbt_(
	    doublereal *, doublereal *, doublereal *, integer *, integer *);
    static integer status;
    extern /* Subroutine */ int fvzeroa_(doublereal *, integer *, integer *);
    extern integer fvcovqq_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), fvluinv_(
	    doublereal *, doublereal *, integer *), fvcovvq_(doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);

/* ----------------------------------------------------------------------- */
/* -- sum up momenta of nt tracks in tList */
/* -- return summed 3-momentum and covariance Dp = cov{p} and Ep = cov{v,p} */
/* -- fill helix parameters and covariance matrix at position new */
/* -- in hl and Ghl */
/* -- */
/* -- restrictions: assumes charge of +/- 1 */
/* --               does not update ql,Cql */
/* -- */
    /* Parameter adjustments */
    --p;
    dp -= 4;
    ep -= 4;
    --tlist;
    --v;
    c__ -= 4;
    ql -= 4;
    cql -= 13;
    hl -= 6;
    ghl -= 31;

    /* Function Body */
    ret_val = 1;
/* -- check for overlap of tList(1..nt) and new */
    i__1 = *nt;
    for (i0 = 1; i0 <= i__1; ++i0) {
	if (*new__ == tlist[i0]) {
	    ret_val = 2;
	}
    }
/* -- zero sum vector and error matrices */
    fvzeroa_(&p[1], &c__3, &c__1);
    fvzeroa_(&dp[4], &c__3, &c__3);
    fvzeroa_(&ep[4], &c__3, &c__3);
/* -- sum up charge, momentum and covariance matrix */
    charge = 0.;
    i__1 = *nt;
    for (i0 = 1; i0 <= i__1; ++i0) {
	it0 = tlist[i0];
	charge += d_sign(&c_b134, &ql[it0 * 3 + 1]);
	fvq2p3_(temp31, temp33, &ql[it0 * 3 + 1], &cql[(it0 * 3 + 1) * 3 + 1],
		 w2pt);
	fvapb_(&p[1], &p[1], temp31, &c__3, &c__1);
	fvapb_(&dp[4], &dp[4], temp33, &c__3, &c__3);
    }
/* -- sum up covariance between different momenta */
    i__1 = *nt;
    for (i0 = 1; i0 <= i__1; ++i0) {
	it0 = tlist[i0];
	status = fvjacobq2p_(j0, &ql[it0 * 3 + 1], w2pt);
	if ((status & 1) != 1) {
	    ret_val = status;
	}
/* -- covariance cov{v,q} */
	status = fvcovvq_(e, &v[1], &c__[4], &ql[it0 * 3 + 1], &ghl[(it0 * 5 
		+ 1) * 5 + 1]);
	if ((status & 1) != 1) {
	    ret_val = status;
	}
	fvabt_(temp33, e, j0, &c__3, &c__3, &c__3);
	fvapb_(&ep[4], &ep[4], temp33, &c__3, &c__3);
/* -- get all Dij = cov{qi,qj}, i != j */
	i__2 = *nt;
	for (i1 = i0 + 1; i1 <= i__2; ++i1) {
	    it1 = tlist[i1];
/* -- sum up correllation terms Dij = cov{qi,qj}: */
/* -- Dp = Dp + sum{i,j,I!=J}( Ji.Dij.Jj^T + Jj.Dij^T.Ji^T) */
	    status = fvcovqq_(d__, &v[1], &c__[4], &ql[it0 * 3 + 1], &ghl[(
		    it0 * 5 + 1) * 5 + 1], &ql[it1 * 3 + 1], &ghl[(it1 * 5 + 
		    1) * 5 + 1]);
	    if ((status & 1) != 1) {
		ret_val = status;
	    }
	    status = fvjacobq2p_(j1, &ql[it1 * 3 + 1], w2pt);
	    if ((status & 1) != 1) {
		ret_val = status;
	    }
	    fvabct_(temp33, j0, d__, j1, &c__3, &c__3);
	    fvapb_(&dp[4], &dp[4], temp33, &c__3, &c__3);
	    fvapbt_(&dp[4], &dp[4], temp33, &c__3, &c__3);
	}
    }
/* -- should fill ql(1,new) and calculate */
/* -- D = dp/dq Dp dp/dq and fill into Cql(1,1,new) */
/* -- but we don't at the moment... */
    fvp2q_(q, &p[1], &charge, w2pt);
    status = fvhch_(&hl[*new__ * 5 + 1], ch, &v[1], q, &c__[4], &dp[4], &ep[4]
	    );
    if ((status & 1) != 1) {
	ret_val = status;
    }
    status = fvluinv_(&ghl[(*new__ * 5 + 1) * 5 + 1], ch, &c__5);
    if ((status & 1) != 1) {
	ret_val = status;
    }
/* CCCC      call fvCopy(Ghl(1,1,new), Ch,dh,dh) */
    return ret_val;
} /* fvsumq_ */

integer fvjacobq2p_(doublereal *j, doublereal *q, doublereal *w2pt)
{
    /* System generated locals */
    integer ret_val;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static doublereal w, cp, tl, sp, pt, psi, dptdw;

/* ----------------------------------------------------------------------- */
/* -- calculate Jacobian for q={w,tl,psi} -> p={px,py,pz} */
/* -- i.e. J(P/Q) */
/* -- */
    /* Parameter adjustments */
    --q;
    j -= 4;

    /* Function Body */
    ret_val = 1;
    w = q[1];
    tl = q[2];
    psi = q[3];
    cp = cos(psi);
    sp = sin(psi);
    pt = *w2pt / abs(w);
    dptdw = -pt / w;
    j[4] = dptdw * cp;
    j[5] = dptdw * sp;
    j[6] = dptdw * tl;
    j[7] = 0.;
    j[8] = 0.;
    j[9] = pt;
    j[10] = -pt * sp;
    j[11] = pt * cp;
    j[12] = 0.;
    return ret_val;
} /* fvjacobq2p_ */

integer fvcovqq_(doublereal *dij, doublereal *v, doublereal *c__, doublereal *
	qi, doublereal *gpi, doublereal *qj, doublereal *gpj)
{
    /* System generated locals */
    integer ret_val;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static doublereal ai[15]	/* was [5][3] */, bi[15]	/* was [5][3] 
	    */, aj[15]	/* was [5][3] */, bj[15]	/* was [5][3] */, wi[
	    9]	/* was [3][3] */, wj[9]	/* was [3][3] */, h0i[5], h0j[5], 
	    temp33[9]	/* was [3][3] */, temp55[25]	/* was [5][5] */;
    extern integer fvabh0_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    extern /* Subroutine */ int fvabct_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *), fvatbc_();
    static doublereal tempp55[25]	/* was [5][5] */;
    static integer status;
    extern integer fvcalcw_(doublereal *, doublereal *, doublereal *);

    /* Fortran I/O blocks */
    static cilist io___81 = { 0, 0, 0, "(1x,a)", 0 };
    static cilist io___83 = { 0, 0, 0, "(1x,a)", 0 };


/* ----------------------------------------------------------------------- */
/* -- calculate covariance matrix Dij between qi and qj */
/* -- */
    /* Parameter adjustments */
    gpj -= 6;
    --qj;
    gpi -= 6;
    --qi;
    c__ -= 4;
    --v;

    /* Function Body */
    ret_val = 1;
/* -- calculate coefficients for measurement equation, A, B, h0 */
    status = fvabh0_(ai, bi, h0i, &v[1], &qi[1]);
    if ((status & 1) != 1) {
	ret_val = status;
    }
    status = fvabh0_(aj, bj, h0j, &v[1], &qj[1]);
    if ((status & 1) != 1) {
	ret_val = status;
    }
/* -- W */
    status = fvcalcw_(wi, &gpi[6], bi);
    if (status == 3) {
	io___81.ciunit = fvcom_1.plun;
	s_wsfe(&io___81);
	do_fio(&c__1, "editing of singular values done during inv. of W", (
		ftnlen)48);
	e_wsfe();
    }
    if ((status & 1) != 1) {
	ret_val = status;
	return ret_val;
    }
    status = fvcalcw_(wj, &gpj[6], bj);
    if (status == 3) {
	io___83.ciunit = fvcom_1.plun;
	s_wsfe(&io___83);
	do_fio(&c__1, "editing of singular values done during inv. of W", (
		ftnlen)48);
	e_wsfe();
    }
    if ((status & 1) != 1) {
	ret_val = status;
	return ret_val;
    }
    fvabct_(temp55, ai, &c__[4], aj, &c__5, &c__3);
    fvabct_(tempp55, &gpi[6], temp55, &gpj[6], &c__5, &c__5);
    fvatbc_(temp33, bi, tempp55, bj, &c__5, &c__3, &c__5);
    fvabct_(dij, wi, temp33, wj, &c__3, &c__3);
    return ret_val;
} /* fvcovqq_ */

integer fvcovvq_(doublereal *e, doublereal *v, doublereal *c__, doublereal *q,
	 doublereal *gh)
{
    /* System generated locals */
    integer ret_val;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static doublereal a[15]	/* was [5][3] */, b[15]	/* was [5][3] */, w[9]
	    	/* was [3][3] */, h0[5], temp33[9]	/* was [3][3] */;
    extern integer fvabh0_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    extern /* Subroutine */ int fvatbc_(), fvabct_(doublereal *, doublereal *,
	     doublereal *, doublereal *, integer *, integer *), fvnega_(
	    doublereal *, doublereal *, integer *, integer *);
    static integer status;
    extern integer fvcalcw_(doublereal *, doublereal *, doublereal *);

    /* Fortran I/O blocks */
    static cilist io___92 = { 0, 0, 0, "(1x,a)", 0 };


/* ----------------------------------------------------------------------- */
/* -- calculate covariance matrix E between v and q */
/* -- */
    /* Parameter adjustments */
    gh -= 6;
    --q;
    c__ -= 4;
    --v;
    e -= 4;

    /* Function Body */
    ret_val = 1;
/* -- calculate coefficients for measurement equation, A, B, h0 */
    status = fvabh0_(a, b, h0, &v[1], &q[1]);
    if ((status & 1) != 1) {
	ret_val = status;
    }
/* -- W */
    status = fvcalcw_(w, &gh[6], b);
    if (status == 3) {
	io___92.ciunit = fvcom_1.plun;
	s_wsfe(&io___92);
	do_fio(&c__1, "editing of singular values done during inv. of W", (
		ftnlen)48);
	e_wsfe();
	ret_val = status;
    }
    if ((status & 1) != 1) {
	ret_val = status;
	return ret_val;
    }
    fvatbc_(temp33, a, &gh[6], b, &c__5, &c__3, &c__5);
    fvabct_(&e[4], &c__[4], temp33, w, &c__3, &c__3);
    fvnega_(&e[4], &e[4], &c__3, &c__3);
    return ret_val;
} /* fvcovvq_ */

integer fvfilter_(doublereal *v, doublereal *c__, doublereal *gv, doublereal *
	q, doublereal *d__, doublereal *e, doublereal *chi2, doublereal *v0, 
	doublereal *gv0, doublereal *p0, doublereal *gp0)
{
    /* Initialized data */

    static integer maxfiter = 10;
    static logical flag__ = TRUE_;

    /* System generated locals */
    integer ret_val;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    extern integer fvsvdfit_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    static doublereal a[15]	/* was [5][3] */, b[15]	/* was [5][3] */, h0[
	    5], q00[3], v00[3];
    extern integer fvfilterer_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *), fvq_(doublereal *, 
	    doublereal *, doublereal *);
    static integer iter;
    extern integer fvabh0_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    extern /* Subroutine */ int fvcopy_(doublereal *, doublereal *, integer *,
	     integer *);
    static integer status;
    static doublereal chi2old;

    /* Fortran I/O blocks */
    static cilist io___104 = { 0, 0, 0, "(1x,a,i5,a,g10.3)", 0 };


/* ----------------------------------------------------------------------- */
/* -- do the Kalman filter step for one track */
/* -- make use of fvFilterer */

/* -- iterate until chi2 change between iterations smaller chi2Cut */
/* -- max #iteration in filter/smoother step */
    /* Parameter adjustments */
    gp0 -= 6;
    --p0;
    gv0 -= 4;
    --v0;
    e -= 4;
    d__ -= 4;
    --q;
    gv -= 4;
    c__ -= 4;
    --v;

    /* Function Body */
    ret_val = 1;
    fvcopy_(v00, &v0[1], &c__3, &c__1);
    status = fvq_(q00, v00, &p0[1]);
    *chi2 = 1e10;
    iter = 0;
L1:
    chi2old = *chi2;
/* -- calculate coefficients for measurement equation, A, B, h0 */
    status = fvabh0_(a, b, h0, v00, q00);
    if ((status & 1) != 1) {
	ret_val = status;
    }
    if (flag__) {
	status = fvfilterer_(&v[1], &c__[4], &gv[4], &q[1], &d__[4], &e[4], 
		chi2, &v0[1], &gv0[4], &p0[1], &gp0[6], a, b, h0);
    } else {
	status = fvsvdfit_(&v[1], &c__[4], &gv[4], &q[1], &d__[4], &e[4], 
		chi2, &v0[1], &gv0[4], &p0[1], &gp0[6], a, b, h0);
    }
    if ((status & 1) != 1) {
	ret_val = status;
    }
/* -- use resultant vertex and 3-mom. for recalculating derivative */
/* -- and do it again */
    fvcopy_(v00, &v[1], &c__3, &c__1);
    fvcopy_(q00, &q[1], &c__3, &c__1);
    ++iter;
/* -- print out chi2 for each iteration */
    if (fvcom_1.print) {
	io___104.ciunit = fvcom_1.plun;
	s_wsfe(&io___104);
	do_fio(&c__1, "Filter iteration ", (ftnlen)17);
	do_fio(&c__1, (char *)&iter, (ftnlen)sizeof(integer));
	do_fio(&c__1, " yields chi2 ", (ftnlen)13);
	do_fio(&c__1, (char *)&(*chi2), (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    if (iter < maxfiter && (d__1 = chi2old - *chi2, abs(d__1)) > .5) {
	goto L1;
    }
    return ret_val;
} /* fvfilter_ */

integer fvsmooth_(doublereal *q, doublereal *d__, doublereal *e, doublereal *
	chi2, doublereal *v, doublereal *c__, doublereal *q0, doublereal *h__,
	 doublereal *gh)
{
    /* Initialized data */

    static integer maxsiter = 1;

    /* System generated locals */
    integer ret_val;
    doublereal d__1;

    /* Local variables */
    static doublereal a[15]	/* was [5][3] */, b[15]	/* was [5][3] */, h0[
	    5];
    extern integer fvsmoother_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    static integer iter;
    extern integer fvabh0_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    extern /* Subroutine */ int fvcopy_(doublereal *, doublereal *, integer *,
	     integer *);
    static integer status;
    static doublereal chi2old;

/* ----------------------------------------------------------------------- */
/* -- calculated smoothed momentum q and cov. matrix D,E */
/* -- D = cov{q} and E = cov{v,q} */
/* -- with smoothed vertex position v,C */
/* -- use filtered momentum q0 = {w,tl,psi} */
/* -- return smoothed chi-square chi2 */

/* -- iterate until chi2 change between iterations smaller chi2Cut */
/* -- iterate only once at the moment */
    /* Parameter adjustments */
    gh -= 6;
    --h__;
    --q0;
    c__ -= 4;
    --v;
    e -= 4;
    d__ -= 4;
    --q;

    /* Function Body */
    ret_val = 1;
/* -- use filtered q0 as expansion point for measurement equation */
    fvcopy_(&q[1], &q0[1], &c__3, &c__1);
    *chi2 = 1e10;
    iter = 0;
L1:
    chi2old = *chi2;
/* -- calculate coefficients for measurement equation, A, B, h0 */
    status = fvabh0_(a, b, h0, &v[1], &q[1]);
    if ((status & 1) != 1) {
	ret_val = status;
    }
    status = fvsmoother_(&q[1], &d__[4], &e[4], chi2, &v[1], &c__[4], &h__[1],
	     &gh[6], a, b, h0);
    if ((status & 1) != 1) {
	ret_val = status;
    }
/* -- check on difference in chi2, ev. re-iterate */
/* -- (maxSIter is 1 for the moment) */
    ++iter;
    if (iter < maxsiter && (d__1 = chi2old - *chi2, abs(d__1)) > .5) {
	goto L1;
    }
    return ret_val;
} /* fvsmooth_ */

integer fvremove_(doublereal *vp, doublereal *cvp, doublereal *gvp, 
	doublereal *ql, doublereal *cql, doublereal *chi2l, integer *nt, 
	integer *tlist, doublereal *alpha, doublereal *v, doublereal *gv, 
	doublereal *hl, doublereal *ghl)
{
    /* System generated locals */
    integer ret_val, i__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    , s_wsfe(cilist *), e_wsfe(void);

    /* Local variables */
    extern integer fvretlif_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), fvsmooth_(doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal a[15]	/* was [5][3] */, b[15]	/* was [5][3] */, e[9]
	    	/* was [3][3] */;
    static integer i__;
    static doublereal h0[5], v0[3];
    static integer it;
    static doublereal gv0[9]	/* was [3][3] */, chi2;
    static integer imax;
    static doublereal prob;
    extern integer fvabh0_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    extern doublereal fvprob_(doublereal *, integer *);
    extern /* Subroutine */ int fvcopy_(doublereal *, doublereal *, integer *,
	     integer *);
    static integer status;
    static doublereal chi2max;
    static logical removed;
    extern /* Subroutine */ int fverror_(char *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___123 = { 0, fvcom_1.text, 0, "(a,i4)", 132, 1 };
    static icilist io___125 = { 0, fvcom_1.text, 0, "(a,i4)", 132, 1 };
    static cilist io___127 = { 0, 0, 0, "(1x,a,i5,a,g10.3,a,g10.3)", 0 };
    static icilist io___128 = { 0, fvcom_1.text, 0, "(a,i4)", 132, 1 };
    static icilist io___130 = { 0, fvcom_1.text, 0, "(a,i4)", 132, 1 };


/* ----------------------------------------------------------------------- */
/* -- check for outliers and remove them from tList */
/* -- */
/* -- calculate chi-square for each track to belong to vertex v */
/* -- remove tracks with chi-square probability larger than alpha */
/* -- return new vertex vp with covariance matrix Cvp, Gvp */
/* -- and updated 3-momenta ql, Cql with chi-square in chi2l */
/* -- removed tracks get a 0 in tList */

    /* Parameter adjustments */
    --vp;
    cvp -= 4;
    gvp -= 4;
    ql -= 4;
    cql -= 13;
    --chi2l;
    --tlist;
    --v;
    gv -= 4;
    hl -= 6;
    ghl -= 31;

    /* Function Body */
    ret_val = 1;
    removed = FALSE_;
/* -- do not do remove with two tracks only */
/* -- set chi**2 to 1 */
    if (*nt <= 2) {
	it = tlist[1];
	if (it != 0) {
	    chi2l[it] = 2.;
	}
	it = tlist[2];
	if (it != 0) {
	    chi2l[it] = 2.;
	}
	fvcopy_(&vp[1], &v[1], &c__3, &c__1);
	fvcopy_(&cvp[4], &gv[4], &c__3, &c__3);
	fvcopy_(&gvp[4], &gv[4], &c__3, &c__3);
	return ret_val;
    }
    fvcopy_(v0, &v[1], &c__3, &c__1);
    fvcopy_(gv0, &gv[4], &c__3, &c__3);
L1:
    imax = 0;
    chi2max = 0.;
    i__1 = *nt;
    for (i__ = 1; i__ <= i__1; ++i__) {
	it = tlist[i__];
/* -- check if this track has been removed */
	if (it != 0) {
/* -- calculate coefficients for measurement equation, A, B, h0 */
/* -- at smoothed vertex position v and smoothed 3-momentum q */
	    status = fvabh0_(a, b, h0, v0, &ql[it * 3 + 1]);
	    if ((status & 1) != 1) {
		s_wsfi(&io___123);
		do_fio(&c__1, "Problems in fvABh0 (outlier) for Track ", (
			ftnlen)39);
		do_fio(&c__1, (char *)&it, (ftnlen)sizeof(integer));
		e_wsfi();
		fverror_(fvcom_1.text, (ftnlen)132);
		ret_val = status;
	    }
/* -- calculate vertex vp, Cvp w/o track p0 and chi2 */
	    status = fvretlif_(&vp[1], &cvp[4], &gvp[4], &chi2, v0, gv0, &ql[
		    it * 3 + 1], &hl[it * 5 + 1], &ghl[(it * 5 + 1) * 5 + 1], 
		    a, b, h0);
	    if ((status & 1) != 1) {
		s_wsfi(&io___125);
		do_fio(&c__1, "Problems in fvRetlif for Track ", (ftnlen)31);
		do_fio(&c__1, (char *)&it, (ftnlen)sizeof(integer));
		e_wsfi();
		fverror_(fvcom_1.text, (ftnlen)132);
		ret_val = status;
	    }
/* -- save this chi2 in chi2l(it) and save i with max. chi2 */
	    chi2l[it] = chi2;
	    if (chi2 > chi2max) {
		chi2max = chi2;
		imax = i__;
	    }
	    if (fvcom_1.print) {
		prob = fvprob_(&chi2, &c__2);
		io___127.ciunit = fvcom_1.plun;
		s_wsfe(&io___127);
		do_fio(&c__1, "Outlier Test for Track ", (ftnlen)23);
		do_fio(&c__1, (char *)&it, (ftnlen)sizeof(integer));
		do_fio(&c__1, " chi2 ", (ftnlen)6);
		do_fio(&c__1, (char *)&chi2, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, " prob. ", (ftnlen)7);
		do_fio(&c__1, (char *)&prob, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	}
/* track has not been removed */
    }
/* -- remove track with largest chi2 if prob(chi2) < alpha */
    if (fvprob_(&chi2max, &c__2) < *alpha) {
	it = tlist[imax];
/* -- remove track iMax from vertex and update all momenta */
	status = fvretlif_(&vp[1], &cvp[4], &gvp[4], &chi2, v0, gv0, &ql[it * 
		3 + 1], &hl[it * 5 + 1], &ghl[(it * 5 + 1) * 5 + 1], a, b, h0)
		;
	if ((status & 1) != 1) {
	    s_wsfi(&io___128);
	    do_fio(&c__1, "Problems in fvRetlif for Track ", (ftnlen)31);
	    do_fio(&c__1, (char *)&it, (ftnlen)sizeof(integer));
	    e_wsfi();
	    fverror_(fvcom_1.text, (ftnlen)132);
	    ret_val = status;
	}
	tlist[imax] = 0;
	i__1 = *nt;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    it = tlist[i__];
	    if (it != 0) {
		status = fvsmooth_(&ql[it * 3 + 1], &cql[(it * 3 + 1) * 3 + 1]
			, e, &chi2, &vp[1], &cvp[4], &ql[it * 3 + 1], &hl[it *
			 5 + 1], &ghl[(it * 5 + 1) * 5 + 1]);
		if ((status & 1) != 1) {
		    s_wsfi(&io___130);
		    do_fio(&c__1, "Problems in fvSmooth for Track ", (ftnlen)
			    31);
		    do_fio(&c__1, (char *)&it, (ftnlen)sizeof(integer));
		    e_wsfi();
		    fverror_(fvcom_1.text, (ftnlen)132);
		    ret_val = status;
		}
	    }
	}
/* - do it again */
	goto L1;
    }
    return ret_val;
} /* fvremove_ */

integer fvsvdfit_(doublereal *v, doublereal *c__, doublereal *gv, doublereal *
	q, doublereal *d__, doublereal *e, doublereal *chi2, doublereal *v0, 
	doublereal *gv0, doublereal *h__, doublereal *gh, doublereal *a, 
	doublereal *b, doublereal *h0)
{
    /* System generated locals */
    integer ret_val;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    extern integer fvsvdinv_(doublereal *, integer *);
    static doublereal m[5], dm[5], bq[3], bv[3];
    extern /* Subroutine */ int fvab_(doublereal *, doublereal *, doublereal *
	    , integer *, integer *, integer *), fvamb_(doublereal *, 
	    doublereal *, doublereal *, integer *, integer *), fvapb_(
	    doublereal *, doublereal *, doublereal *, integer *, integer *), 
	    fvatb_(doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, integer *);
    static doublereal temp31[3], temp33[9]	/* was [3][3] */, temp51[5], 
	    chi2ar[1]	/* was [1][1] */;
    extern /* Subroutine */ int fvatbc_();
    static doublereal tempp31[3];
    extern /* Subroutine */ int fvcopy_(doublereal *, doublereal *, integer *,
	     integer *);
    static integer status;
    extern integer fvinv6s_(doublereal *, doublereal *, doublereal *);
    extern /* Subroutine */ int fvsatba_(doublereal *, doublereal *, 
	    doublereal *, integer *, integer *);

    /* Fortran I/O blocks */
    static cilist io___133 = { 0, 0, 0, "(1x,a)", 0 };
    static cilist io___134 = { 0, 0, 0, "(1x,a)", 0 };


/* -------------------------------------------------------------------------- */
/* -- calculate v, q, C, D, E, chi2 by solving the normal equation */
/* -- */
/* -- / Gv0 + A^T.Gh.A   A^T.Gh.B \  / v \     / Gv0.v0 + A^T.Gh.(h-h0) \ */
/* -- |                            | |    | =  |                         | */
/* -- \    B^T.Gh.A      B^T.Gh.B /  \ q /     \          B^T.Gh.(h-h0) / */
/* -- */
/* -- Ansatz: */
/* -- */
/* -- / v \    / C  E \    / Gv0.v0 + A^T.Gh.(h-h0) \ */
/* -- |    | = |  T    | . |                         | */
/* -- \ q /    \ E  D /    \          B^T.Gh.(h-h0) / */
/* -- */
    /* Parameter adjustments */
    --h0;
    b -= 6;
    a -= 6;
    gh -= 6;
    --h__;
    gv0 -= 4;
    --v0;
    e -= 4;
    d__ -= 4;
    --q;
    gv -= 4;
    c__ -= 4;
    --v;

    /* Function Body */
    ret_val = 1;
/* -- C <- Gv0 + A^T.Gh.A */
    fvsatba_(temp33, &a[6], &gh[6], &c__5, &c__3);
    fvapb_(&c__[4], &gv0[4], temp33, &c__3, &c__3);
/* -- D <- B^T.Gh.B */
    fvsatba_(&d__[4], &b[6], &gh[6], &c__5, &c__3);
/* -- E <- A^T.Gh.B */
    fvatbc_(&e[4], &a[6], &gh[6], &b[6], &c__5, &c__3);
    status = fvinv6s_(&c__[4], &d__[4], &e[4]);
    if (status == 3) {
	io___133.ciunit = fvcom_1.plun;
	s_wsfe(&io___133);
	do_fio(&c__1, "editing of singular values done during inv. of CDE", (
		ftnlen)50);
	e_wsfe();
    }
    if ((status & 1) != 1) {
	ret_val = status;
	return ret_val;
    }
/* -- Gv = C^(-1) */
    fvcopy_(&gv[4], &c__[4], &c__3, &c__3);
    status = fvsvdinv_(&gv[4], &c__3);
    if (status == 3) {
	io___134.ciunit = fvcom_1.plun;
	s_wsfe(&io___134);
	do_fio(&c__1, "editing of singular values done inverting C in fvSVDf"
		"it", (ftnlen)55);
	e_wsfe();
    }
    if ((status & 1) != 1) {
	ret_val = status;
    }
/* -- m = h - h0 */
    fvamb_(m, &h__[1], &h0[1], &c__5, &c__1);
/* -- bv = Gv0.v0 + A^T.Gh.(h - h0) */
/* -- bq = B^T.Gh.(h - h0) */
    fvab_(temp31, &gv0[4], &v0[1], &c__3, &c__3, &c__1);
    fvab_(temp51, &gh[6], m, &c__5, &c__5, &c__1);
    fvatb_(tempp31, &a[6], temp51, &c__5, &c__3, &c__1);
    fvapb_(bv, temp31, tempp31, &c__3, &c__1);
    fvatb_(bq, &b[6], temp51, &c__5, &c__3, &c__1);
/* -- v = C.bv + E.bq */
    fvab_(&v[1], &c__[4], bv, &c__3, &c__3, &c__1);
    fvab_(temp31, &e[4], bq, &c__3, &c__3, &c__1);
    fvapb_(&v[1], &v[1], temp31, &c__3, &c__1);
/* -- q = E^T.bv + D.bq */
    fvatb_(&q[1], &e[4], bv, &c__3, &c__3, &c__1);
    fvab_(temp31, &d__[4], bq, &c__3, &c__3, &c__1);
    fvapb_(&q[1], &q[1], temp31, &c__3, &c__1);
/* -- dm = h - h0 - A.v */
    fvab_(temp51, &a[6], &v[1], &c__5, &c__3, &c__1);
    fvamb_(dm, m, temp51, &c__5, &c__1);
/* -- chi2 = (h - h0 - A.v - B.q)^T.Gh.(h - h0 - A.v - B.q) + */
/* --        (v-v0)^T.Gv0.(v-v0) */
    fvab_(temp51, &b[6], &q[1], &c__5, &c__3, &c__1);
    fvamb_(temp51, dm, temp51, &c__5, &c__1);
    fvsatba_(chi2ar, temp51, &gh[6], &c__5, &c__1);
    *chi2 = chi2ar[0];
    fvamb_(temp31, &v[1], &v0[1], &c__3, &c__1);
    fvsatba_(chi2ar, temp31, &gv0[4], &c__3, &c__1);
    *chi2 += chi2ar[0];
    return ret_val;
} /* fvsvdfit_ */

integer fvfilterer_(doublereal *v, doublereal *c__, doublereal *gv, 
	doublereal *q, doublereal *d__, doublereal *e, doublereal *chi2, 
	doublereal *v0, doublereal *gv0, doublereal *p0, doublereal *g, 
	doublereal *a, doublereal *b, doublereal *h0)
{
    /* System generated locals */
    integer ret_val;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    , s_wsfe(cilist *), e_wsfe(void);

    /* Local variables */
    extern /* Subroutine */ int fvcalcgb_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    extern integer fvsvdinv_(doublereal *, integer *);
    static doublereal m[5], w[9]	/* was [3][3] */, gb[25]	/* 
	    was [5][5] */, dm[5];
    extern /* Subroutine */ int fvab_(doublereal *, doublereal *, doublereal *
	    , integer *, integer *, integer *), fvamb_(doublereal *, 
	    doublereal *, doublereal *, integer *, integer *), fvapb_(
	    doublereal *, doublereal *, doublereal *, integer *, integer *), 
	    fvatb_(doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, integer *);
    static doublereal temp31[3], temp33[9]	/* was [3][3] */, temp51[5], 
	    temp55[25]	/* was [5][5] */, chi2ar[1]	/* was [1][1] */;
    extern /* Subroutine */ int fvatbc_(), fvnega_(doublereal *, doublereal *,
	     integer *, integer *);
    static doublereal tempp31[3], tempp55[25]	/* was [5][5] */;
    extern /* Subroutine */ int fvcopy_(doublereal *, doublereal *, integer *,
	     integer *);
    static integer status;
    extern integer fvcalcw_(doublereal *, doublereal *, doublereal *);
    extern /* Subroutine */ int fvsatba_(doublereal *, doublereal *, 
	    doublereal *, integer *, integer *), fvsabat_(doublereal *, 
	    doublereal *, doublereal *, integer *, integer *), fverror_(char *
	    , ftnlen);

    /* Fortran I/O blocks */
    static icilist io___145 = { 0, fvcom_1.text, 0, "(a)", 132, 1 };
    static cilist io___147 = { 0, 0, 0, "(1x,a)", 0 };


/* -------------------------------------------------------------------------- */
/* -- Kalman filter for vertex fit */
/* -- */
/* -- needs 2 inversions of 3x3-matrices */
/* -- */
/* -- author: Lothar A.T. Bauerdick, CERN/PPE */
/* -- date:   Jan. 1991 */
/* -- */
    /* Parameter adjustments */
    --h0;
    b -= 6;
    a -= 6;
    g -= 6;
    --p0;
    gv0 -= 4;
    --v0;
    e -= 4;
    d__ -= 4;
    --q;
    gv -= 4;
    c__ -= 4;
    --v;

    /* Function Body */
    ret_val = 1;
/* -- W */
    status = fvcalcw_(w, &g[6], &b[6]);
    if (status == 3) {
	s_wsfi(&io___145);
	do_fio(&c__1, "editing of singular values done during inv. of W", (
		ftnlen)48);
	e_wsfi();
	fverror_(fvcom_1.text, (ftnlen)132);
/* -- ???? */
/* -- for now: forget about the fit... */
	ret_val = 2;
	return ret_val;
    }
    if ((status & 1) != 1) {
	ret_val = status;
	return ret_val;
    }
/* -- GB */
    fvcalcgb_(gb, w, &b[6], &g[6]);
/* -- C */
    fvsatba_(&gv[4], &a[6], gb, &c__5, &c__3);
    fvapb_(&gv[4], &gv0[4], &gv[4], &c__3, &c__3);
    fvcopy_(&c__[4], &gv[4], &c__3, &c__3);
    status = fvsvdinv_(&c__[4], &c__3);
    if (status == 3) {
	io___147.ciunit = fvcom_1.plun;
	s_wsfe(&io___147);
	do_fio(&c__1, "editing of singular values done during inv. of C^(-1)",
		 (ftnlen)53);
	e_wsfe();
    }
    if ((status & 1) != 1) {
	ret_val = status;
	return ret_val;
    }
/* -- m = p0 - h0 */
    fvamb_(m, &p0[1], &h0[1], &c__5, &c__1);
/* -- v */
    fvab_(temp51, gb, m, &c__5, &c__5, &c__1);
    fvatb_(temp31, &a[6], temp51, &c__5, &c__3, &c__1);
    fvab_(tempp31, &gv0[4], &v0[1], &c__3, &c__3, &c__1);
    fvapb_(temp31, tempp31, temp31, &c__3, &c__1);
    fvab_(&v[1], &c__[4], temp31, &c__3, &c__3, &c__1);
/* -- dm = p0 - h0 - A.v */
    fvab_(temp51, &a[6], &v[1], &c__5, &c__3, &c__1);
    fvamb_(dm, m, temp51, &c__5, &c__1);
/* -- q */
    fvab_(temp51, &g[6], dm, &c__5, &c__5, &c__1);
    fvatb_(temp31, &b[6], temp51, &c__5, &c__3, &c__1);
    fvab_(&q[1], w, temp31, &c__3, &c__3, &c__1);
/* -- D */
    fvsabat_(temp55, &a[6], &c__[4], &c__5, &c__3);
    fvsatba_(tempp55, &g[6], temp55, &c__5, &c__5);
    fvsatba_(temp33, &b[6], tempp55, &c__5, &c__3);
    fvsatba_(&d__[4], w, temp33, &c__3, &c__3);
    fvapb_(&d__[4], w, &d__[4], &c__3, &c__3);
/* -- E */
    fvatbc_(temp33, &b[6], &g[6], &a[6], &c__5, &c__3);
    fvatbc_(&e[4], w, temp33, &c__[4], &c__3, &c__3);
    fvnega_(&e[4], &e[4], &c__3, &c__3);
/* -- chi2 */
    fvab_(temp51, &b[6], &q[1], &c__5, &c__3, &c__1);
    fvamb_(temp51, dm, temp51, &c__5, &c__1);
    fvsatba_(chi2ar, temp51, &g[6], &c__5, &c__1);
    *chi2 = chi2ar[0];
    fvamb_(temp31, &v[1], &v0[1], &c__3, &c__1);
    fvsatba_(chi2ar, temp31, &gv0[4], &c__3, &c__1);
    *chi2 += chi2ar[0];
    return ret_val;
} /* fvfilterer_ */

integer fvsmoother_(doublereal *q, doublereal *d__, doublereal *e, doublereal 
	*chi2, doublereal *v, doublereal *c__, doublereal *h__, doublereal *g,
	 doublereal *a, doublereal *b, doublereal *h0)
{
    /* System generated locals */
    integer ret_val;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    extern /* Subroutine */ int fvcalcgb_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal m[5], w[9]	/* was [3][3] */, gb[25]	/* 
	    was [5][5] */, dm[5];
    extern integer fvh_(doublereal *, doublereal *, doublereal *);
    extern /* Subroutine */ int fvab_(doublereal *, doublereal *, doublereal *
	    , integer *, integer *, integer *), fvamb_(doublereal *, 
	    doublereal *, doublereal *, integer *, integer *), fvapb_(
	    doublereal *, doublereal *, doublereal *, integer *, integer *), 
	    fvatb_(doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, integer *);
    static doublereal temp31[3], temp33[9]	/* was [3][3] */, temp51[5], 
	    temp55[25]	/* was [5][5] */, chi2ar[1]	/* was [1][1] */;
    extern /* Subroutine */ int fvatbc_(), fvnega_(doublereal *, doublereal *,
	     integer *, integer *);
    static doublereal tempp55[25]	/* was [5][5] */;
    static integer status;
    extern integer fvcalcw_(doublereal *, doublereal *, doublereal *);
    extern /* Subroutine */ int fvsabat_(doublereal *, doublereal *, 
	    doublereal *, integer *, integer *), fvsatba_(doublereal *, 
	    doublereal *, doublereal *, integer *, integer *);

    /* Fortran I/O blocks */
    static cilist io___159 = { 0, 0, 0, "(1x,a)", 0 };


/* -------------------------------------------------------------------------- */
/* -- smoother */
/* -- */
/* -- author: Lothar A.T. Bauerdick, CERN/PPE */
/* -- date:   Jan. 1991 */
/* -- */
    /* Parameter adjustments */
    --h0;
    b -= 6;
    a -= 6;
    g -= 6;
    --h__;
    c__ -= 4;
    --v;
    e -= 4;
    d__ -= 4;
    --q;

    /* Function Body */
    ret_val = 1;
/* -- W */
    status = fvcalcw_(w, &g[6], &b[6]);
    if (status == 3) {
	io___159.ciunit = fvcom_1.plun;
	s_wsfe(&io___159);
	do_fio(&c__1, "editing of singular values done during inv. of W", (
		ftnlen)48);
	e_wsfe();
    }
    if ((status & 1) != 1) {
	ret_val = status;
	return ret_val;
    }
/* -- GB */
    fvcalcgb_(gb, w, &b[6], &g[6]);
/* -- m = h - h0 */
    fvamb_(m, &h__[1], &h0[1], &c__5, &c__1);
/* -- dm = h - h0 - A.v */
    fvab_(temp51, &a[6], &v[1], &c__5, &c__3, &c__1);
    fvamb_(dm, m, temp51, &c__5, &c__1);
/* -- q */
    fvab_(temp51, &g[6], dm, &c__5, &c__5, &c__1);
    fvatb_(temp31, &b[6], temp51, &c__5, &c__3, &c__1);
    fvab_(&q[1], w, temp31, &c__3, &c__3, &c__1);
/* -- D */
    fvsabat_(temp55, &a[6], &c__[4], &c__5, &c__3);
    fvsatba_(tempp55, &g[6], temp55, &c__5, &c__5);
    fvsatba_(temp33, &b[6], tempp55, &c__5, &c__3);
    fvsatba_(&d__[4], w, temp33, &c__3, &c__3);
    fvapb_(&d__[4], w, &d__[4], &c__3, &c__3);
/* -- E */
    fvatbc_(temp33, &b[6], &g[6], &a[6], &c__5, &c__3);
    fvatbc_(&e[4], w, temp33, &c__[4], &c__3, &c__3);
    fvnega_(&e[4], &e[4], &c__3, &c__3);
/* -- chi2 = (h - h(v,q))^T.Gh(v,q).(h - h(v,q)) */
/* -- where we calculate Gh(v,q) from the fit-result covariance */
/* -- matrices C,D,E for the smoothed results v,q */
/* CCCC      call fvAB(temp51, B,q,dh,dv,1) */
/* CCCC      call fvAMB(temp51, dm,temp51,dh,1) */
/* CC      status = fvCh(Gh, v,q,C,D,E) */
/* CC      if (status .NE. NORMAL) then */
/* CC         fvSmoother = status */
/* CC      endif */
/* CC      status = fvSVDinv(Gh,dh) */
/* CC      if (status .NE. NORMAL) then */
/* CC         fvSmoother = status */
/* CC      endif */
    status = fvh_(temp51, &v[1], &q[1]);
    if ((status & 1) != 1) {
	ret_val = status;
    }
    fvamb_(temp51, &h__[1], temp51, &c__5, &c__1);
/* CC      call fvsATBA(chi2ar, temp51,Gh,dh,1) */
    fvsatba_(chi2ar, temp51, &g[6], &c__5, &c__1);
    *chi2 = chi2ar[0];
    return ret_val;
} /* fvsmoother_ */

integer fvretlif_(doublereal *vp, doublereal *cp, doublereal *gvp, doublereal 
	*chi2, doublereal *v, doublereal *gv, doublereal *q, doublereal *p0, 
	doublereal *g, doublereal *a, doublereal *b, doublereal *h0)
{
    /* System generated locals */
    integer ret_val;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    extern /* Subroutine */ int fvcalcgb_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    extern integer fvsvdinv_(doublereal *, integer *);
    static doublereal m[5], w[9]	/* was [3][3] */, gb[25]	/* 
	    was [5][5] */, dm[5];
    extern /* Subroutine */ int fvab_(doublereal *, doublereal *, doublereal *
	    , integer *, integer *, integer *), fvamb_(doublereal *, 
	    doublereal *, doublereal *, integer *, integer *), fvatb_(
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *);
    static doublereal temp31[3], temp51[5], chi2ar[1]	/* was [1][1] */, 
	    tempp31[3];
    extern /* Subroutine */ int fvcopy_(doublereal *, doublereal *, integer *,
	     integer *);
    static integer status;
    extern integer fvcalcw_(doublereal *, doublereal *, doublereal *);
    extern /* Subroutine */ int fvsatba_(doublereal *, doublereal *, 
	    doublereal *, integer *, integer *);

    /* Fortran I/O blocks */
    static cilist io___171 = { 0, 0, 0, "(1x,a)", 0 };
    static cilist io___173 = { 0, 0, 0, "(1x,a)", 0 };


/* -------------------------------------------------------------------------- */
/* -- inverse Kalman filter, removes track p0, G from vertex v, Gv */
/* -- returns new vertex vp, Cp and the chi2 for track belonging to */
/* -- new vertex */
/* -- */
/* -- author: Lothar A.T. Bauerdick, CERN/PPE */
/* -- date:   Jan. 1991 */
/* -- */
    /* Parameter adjustments */
    --h0;
    b -= 6;
    a -= 6;
    g -= 6;
    --p0;
    --q;
    gv -= 4;
    --v;
    gvp -= 4;
    cp -= 4;
    --vp;

    /* Function Body */
    ret_val = 1;
/* -- W */
    status = fvcalcw_(w, &g[6], &b[6]);
    if (status == 3) {
	io___171.ciunit = fvcom_1.plun;
	s_wsfe(&io___171);
	do_fio(&c__1, "editing of singular values done during inv. of W", (
		ftnlen)48);
	e_wsfe();
    }
    if ((status & 1) != 1) {
	ret_val = status;
	return ret_val;
    }
/* -- GB */
    fvcalcgb_(gb, w, &b[6], &g[6]);
/* -- Gvp */
    fvsatba_(&gvp[4], &a[6], gb, &c__5, &c__3);
    fvamb_(&gvp[4], &gv[4], &gvp[4], &c__3, &c__3);
/* -- Cp = Gvp^(-1) */
    fvcopy_(&cp[4], &gvp[4], &c__3, &c__3);
    status = fvsvdinv_(&cp[4], &c__3);
    if (status == 3) {
	io___173.ciunit = fvcom_1.plun;
	s_wsfe(&io___173);
	do_fio(&c__1, "editing of singular values done during inv. of C^(-1)",
		 (ftnlen)53);
	e_wsfe();
    }
    if ((status & 1) != 1) {
	ret_val = status;
	return ret_val;
    }
/* -- m = p0 - h0 */
    fvamb_(m, &p0[1], &h0[1], &c__5, &c__1);
/* -- vp = Cp.(Gv.v - A^T.GB.(p0-h0)) */
    fvab_(temp31, &gv[4], &v[1], &c__3, &c__3, &c__1);
    fvab_(temp51, gb, m, &c__5, &c__5, &c__1);
    fvatb_(tempp31, &a[6], temp51, &c__5, &c__3, &c__1);
    fvamb_(temp31, temp31, tempp31, &c__3, &c__1);
    fvab_(&vp[1], &cp[4], temp31, &c__3, &c__3, &c__1);
/* -- chi2 */

/* -- chi2 = r^T.G.r + (v-vp)^T.Gvp.(v-vp) */
/* -- r = p0 - h0 - A.v - B.q */
/* -- BUT: */
/* -- shoudn't this be h(vp,qp) at new vertex vp? */
/* -- i.e. should q be calculated using fvq(q, vp,p0)? */
/* -- do it... */
/* CC      call fvq(qp, vp,p0) */
/* -- Fruehwirt CERN 90-06 says: */
/* -- distance of track from the new vertex is expressed by the chi-square */
/* -- of the smoothed residuals */
/* -- chi2 = (p0-h0-A.v-B.q)^T.G.(p0-h0-A.v-B.q) + (v-vp)^T.Gvp.(v-vp) */
    fvab_(temp51, &a[6], &v[1], &c__5, &c__3, &c__1);
    fvamb_(dm, m, temp51, &c__5, &c__1);
    fvab_(temp51, &b[6], &q[1], &c__5, &c__3, &c__1);
    fvamb_(temp51, dm, temp51, &c__5, &c__1);
    fvsatba_(chi2ar, temp51, &g[6], &c__5, &c__1);
    *chi2 = chi2ar[0];
    fvamb_(temp31, &v[1], &vp[1], &c__3, &c__1);
    fvsatba_(chi2ar, temp31, &gvp[4], &c__3, &c__1);
    *chi2 += chi2ar[0];
    return ret_val;
} /* fvretlif_ */

integer fvv0q0_(doublereal *v0, doublereal *cv0, doublereal *q0, doublereal *
	d0, doublereal *x0, doublereal *cx0, doublereal *p0, doublereal *cp0)
{
    /* System generated locals */
    integer ret_val;

    /* Local variables */
    extern integer fvq_(doublereal *, doublereal *, doublereal *), fvvr_(
	    doublereal *, doublereal *), fvcvr_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static integer status;

/* -------------------------------------------------------------------------- */
/* -- calculate and return v0, Cv0, q0, D0 from x0, Cx0, p0, Cp0 */
/* -- */
/* -- calculate vertex position v0 and 3-momentum vector q0 */
/* -- from track parameters p0 (reference surface: r = d0) */
/* -- at carthesian coordinates x0 */
/* -- */
/* -- p0 = {w, tl, psi0, d0, z0}, x0 = {x, y, z} */
/* -- */
/* -- ???? */
/* -- we do not actually calculate D0, as it is not used */
/* -- otherwise it should be 1/delta, delta -> inf. */
/* -- ???? */
/* -- */

    /* Parameter adjustments */
    cp0 -= 6;
    --p0;
    cx0 -= 4;
    --x0;
    d0 -= 4;
    --q0;
    cv0 -= 4;
    --v0;

    /* Function Body */
    ret_val = 1;
/* -- calculate vertex vector in appropriate coord. system */
    status = fvvr_(&v0[1], &x0[1]);
    if ((status & 1) != 1) {
	ret_val = status;
	return ret_val;
    }
/* -- calculate direction vector q0 at min dist. of helix p0 to point v0 */
    status = fvq_(&q0[1], &v0[1], &p0[1]);
    if ((status & 1) != 1) {
	ret_val = status;
	return ret_val;
    }
/* -- tranform covariance matrix */
    status = fvcvr_(&cv0[4], &v0[1], &x0[1], &cx0[4]);
    if ((status & 1) != 1) {
	ret_val = status;
	return ret_val;
    }
    return ret_val;
} /* fvv0q0_ */

integer fvcvr_(doublereal *cv, doublereal *v, doublereal *x, doublereal *cx)
{
    /* System generated locals */
    integer ret_val;

    /* Local variables */
    static doublereal j[9]	/* was [3][3] */, r__, z__, phi;
    extern /* Subroutine */ int fvsabat_(doublereal *, doublereal *, 
	    doublereal *, integer *, integer *);

/* ---------------------------------------------------------------------- */
/* -- transform covariance matrix Cx of x to get Cv */
/* -- */
    /* Parameter adjustments */
    cx -= 4;
    --x;
    --v;
    cv -= 4;

    /* Function Body */
    ret_val = 1;
    r__ = v[1];
    phi = v[2];
    z__ = v[3];
/* -- Jacobi */
/* -- J = d(r,phi,z)/d(x,y,z) */
    j[0] = (x[1] - 0.) / r__;
    j[3] = x[2] / r__;
    j[6] = 0.;
    j[1] = -x[2] / r__ / r__;
    j[4] = (x[1] - 0.) / r__ / r__;
    j[7] = 0.;
    j[2] = 0.;
    j[5] = 0.;
    j[8] = 1.;
/* -- Cv = J.Cx.J^T */
    fvsabat_(&cv[4], j, &cx[4], &c__3, &c__3);
    return ret_val;
} /* fvcvr_ */

integer fvvr_(doublereal *v, doublereal *x)
{
    /* System generated locals */
    integer ret_val;

    /* Builtin functions */
    double sqrt(doublereal), atan2(doublereal, doublereal);

    /* Local variables */
    static doublereal xx, yy;

/* -------------------------------------------------------------------------- */
/* -- calculate vertex vector v = {r,phi,z} from x = {x, y, z} */
/* -- i.e. transform to a cylindrical coordinate system */
/* -- */

    /* Parameter adjustments */
    --x;
    --v;

    /* Function Body */
    ret_val = 1;
    xx = x[1] - 0.;
    yy = x[2];
    v[1] = sqrt(xx * xx + yy * yy);
    if (v[1] == 0.) {
/* -- can't handle vertex at (0,0,z) at the moment! */
	ret_val = 2;
	return ret_val;
    }
    v[2] = atan2(yy, xx);
    v[3] = x[3];
    return ret_val;
} /* fvvr_ */

integer fvx_(doublereal *x, doublereal *v)
{
    /* System generated locals */
    integer ret_val;

    /* Builtin functions */
    double sin(doublereal), cos(doublereal);

    /* Local variables */
    static doublereal r__, z__, phi, cphi, sphi;

/* -------------------------------------------------------------------------- */
/* -- calculate vector x = {x, y, z} from v = {r,phi,z} */
/* -- i.e. transform back from cylindrical coordinate system */
/* -- */

    /* Parameter adjustments */
    --v;
    --x;

    /* Function Body */
    ret_val = 1;
    r__ = v[1];
    phi = v[2];
    z__ = v[3];
    sphi = sin(phi);
    cphi = cos(phi);
    x[1] = r__ * cphi + 0.;
    x[2] = r__ * sphi;
    x[3] = z__;
    return ret_val;
} /* fvx_ */

integer fvxp_(doublereal *x, doublereal *cx, doublereal *p, doublereal *cp, 
	doublereal *v, doublereal *q, doublereal *c__, doublereal *d__, 
	doublereal *e, doublereal *a, doublereal *b)
{
    /* System generated locals */
    integer ret_val;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal), atan(doublereal);

    /* Local variables */
    static doublereal j[9]	/* was [3][3] */, r__, w, z__, tl, xi, phi, 
	    cxi, psi, sxi, oow, cphi, sphi, gamma;
    extern /* Subroutine */ int fvapb_(doublereal *, doublereal *, doublereal 
	    *, integer *, integer *);
    static doublereal temp55[25]	/* was [5][5] */;
    extern /* Subroutine */ int fvabct_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *), fvsabat_(
	    doublereal *, doublereal *, doublereal *, integer *, integer *), 
	    fvabtct_(doublereal *, doublereal *, doublereal *, doublereal *, 
	    integer *, integer *);

/* -------------------------------------------------------------------------- */
/* -- calculate and return x, Cx, p, Cp from v, q, C, D, E, A, B */
/* -- A and B given for conveniance only */

    /* Parameter adjustments */
    b -= 6;
    a -= 6;
    e -= 4;
    d__ -= 4;
    c__ -= 4;
    --q;
    --v;
    cp -= 6;
    --p;
    cx -= 4;
    --x;

    /* Function Body */
    ret_val = 1;
    r__ = v[1];
    phi = v[2];
    z__ = v[3];
    w = q[1];
    tl = q[2];
    psi = q[3];
    oow = 1. / w;
    xi = psi - phi;
    cxi = cos(xi);
    sxi = sin(xi);
    gamma = atan(r__ * cxi / (oow - r__ * sxi));
    sphi = sin(phi);
    cphi = cos(phi);
    x[1] = r__ * cphi + 0.;
    x[2] = r__ * sphi;
    x[3] = z__;
    p[1] = w;
    p[2] = tl;
    p[3] = psi - gamma;
    p[4] = -r__ * cxi / sin(gamma) + oow;
    p[5] = z__ - gamma / w * tl;
/* -- transform covariance matrix */
/* -- J = d(x,y,z)/d(r,phi,z) */
    j[0] = cphi;
    j[3] = -r__ * sphi;
    j[6] = 0.;
    j[1] = sphi;
    j[4] = r__ * cphi;
    j[7] = 0.;
    j[2] = 0.;
    j[5] = 0.;
    j[8] = 1.;
/* -- Cx = J.C.J^T */
    fvsabat_(&cx[4], j, &c__[4], &c__3, &c__3);
/* -- Cp = A.C.A^T + B.E.A^T + A.E^T.B^T + B.D.B^T */
    fvsabat_(&cp[6], &a[6], &c__[4], &c__5, &c__3);
    fvabct_(temp55, &b[6], &e[4], &a[6], &c__5, &c__3);
    fvapb_(&cp[6], &cp[6], temp55, &c__5, &c__5);
    fvabtct_(temp55, &a[6], &e[4], &b[6], &c__5, &c__3);
    fvapb_(&cp[6], &cp[6], temp55, &c__5, &c__5);
    fvsabat_(temp55, &b[6], &d__[4], &c__5, &c__3);
    fvapb_(&cp[6], &cp[6], temp55, &c__5, &c__5);
    return ret_val;
} /* fvxp_ */

integer fvq_(doublereal *q, doublereal *v, doublereal *h__)
{
    /* System generated locals */
    integer ret_val;

    /* Local variables */
    static doublereal w, tl, psi, psi0, gamma;

/* -------------------------------------------------------------------------- */
/* -- calculate 3-momentum vector q = {w, tl, psi} */
/* -- from helix parameters h (reference surface: r = d0) */
/* -- at minimum distance to space point v */
/* -- */
/* -- h = {w, tl, psi0, d0, z0}, v = {x, y, z} */
/* -- */

/* CC      double precision dd0, r0, u0, xi, gamma, psi, r, phi */
    /* Parameter adjustments */
    --h__;
    --v;
    --q;

    /* Function Body */
    ret_val = 1;
/* CC      r = sqrt(v(1)*v(1) + v(2)*v(2)) */
/* CC      phi =  atan2(v(2),v(1)) */
    w = h__[1];
    tl = h__[2];
    psi0 = h__[3];
/* CC      dd0 =   h(4) */

/* -- calculate psi at a position on track helix p, where */
/* -- the distance between point v and the helix is minimal */
/* -- in the x-y plane. */
/* -- How can we include the distance in z? Is this necessary? */

/* -- for the moment, just return track parameters in global system, */
/* -- i.e. at minimum distance of approach to origin */
/* CC      if (dd0 .EQ. 0d0) then */
    gamma = 0.;
/* CC      else */
/* CC        r0 = abs(dd0) */
/* CC        u0 = 2d0*r0*atan(-dd0*cos(psi0)/(r0+dd0*sin(psi0))) */
/* CC        xi = -u0/dd0 - phi */
/* CC        gamma = atan(r*sin(xi)/(1d0/w-dd0 - r*cos(xi))) */
/* CC      end do */
    psi = psi0 + gamma;
    q[1] = w;
    q[2] = tl;
    q[3] = psi;
    return ret_val;
} /* fvq_ */

integer fvrh2q_(doublereal *q, doublereal *cq, doublereal *r__, doublereal *
	h__, doublereal *ch)
{
    /* System generated locals */
    integer ret_val;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), asin(doublereal);

    /* Local variables */
    static doublereal j[15]	/* was [3][5] */, w, d0;
    static integer i0, i1;
    static doublereal tl, psi, psi0, srmd, sdtw, gamma, dgamdu;
    extern /* Subroutine */ int fvsabat_(doublereal *, doublereal *, 
	    doublereal *, integer *, integer *);

/* -------------------------------------------------------------------------- */
/* -- calculate 3-momentum vector q = {w, tl, psi} and covariance matrix Cq */
/* -- from track parameters h = {w, tl, psi0, d0, z0} (reference surface is */
/* -- Rref = d0) and covariance matrix Ch at a given radius r > d0 */
/* -- */

    /* Parameter adjustments */
    ch -= 6;
    --h__;
    cq -= 4;
    --q;

    /* Function Body */
    ret_val = 1;
    w = h__[1];
    tl = h__[2];
    psi0 = h__[3];
    d0 = h__[4];
/* -- track crosses r only, if */
/* -- d0<r and d0*w < 1 */
    if (abs(d0) < *r__) {
	srmd = sqrt(*r__ * *r__ - d0 * d0);
	sdtw = sqrt(1 - d0 * w);
	gamma = asin(w / 2. * srmd / sdtw) * 2.;
/* Computing 3rd power */
	d__1 = sdtw;
/* Computing 2nd power */
	d__2 = srmd * w / 2. / sdtw;
	dgamdu = .5 / srmd / (d__1 * (d__1 * d__1)) / sqrt(1. - d__2 * d__2);
	psi = psi0 + gamma;
	q[1] = w;
	q[2] = tl;
	q[3] = psi;
	j[0] = 1.;
	j[3] = 0.;
	j[6] = 0.;
	j[9] = 0.;
	j[12] = 0.;
	j[1] = 0.;
	j[4] = 1.;
	j[7] = 0.;
	j[10] = 0.;
	j[13] = 0.;
	j[2] = (2. - d0 * w) * srmd * srmd * dgamdu;
	j[5] = 0.;
	j[8] = 1.;
	j[11] = w * (d0 * -2. + d0 * d0 * w + *r__ * *r__ * w) * dgamdu;
	j[14] = 1.;
	fvsabat_(&cq[4], j, &ch[6], &c__3, &c__5);
    } else {
	psi = psi0;
	for (i0 = 1; i0 <= 3; ++i0) {
	    for (i1 = 1; i1 <= 3; ++i1) {
		cq[i0 + i1 * 3] = ch[i0 + i1 * 5];
	    }
	}
    }
    return ret_val;
} /* fvrh2q_ */

integer fvrh2v_(doublereal *v, doublereal *cv, doublereal *r__, doublereal *
	h__, doublereal *ch)
{
    /* System generated locals */
    integer ret_val;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), asin(doublereal);

    /* Local variables */
    static doublereal j[15]	/* was [3][5] */, w, d0;
    static integer i0, i1;
    static doublereal z0, tl, xx, yy, zz, psi, psi0, srmd, sdtw, gamma, 
	    dgamdu;
    extern /* Subroutine */ int fvsabat_(doublereal *, doublereal *, 
	    doublereal *, integer *, integer *);

/* -------------------------------------------------------------------------- */
/* -- calculate vertex vector v = {vx, vy, vz} and covariance matrix Cv */
/* -- from track parameters h = {w, tl, psi0, d0, z0} (reference surface is */
/* -- Rref = d0) and covariance matrix Ch at a given radius r > d0 */

/* ---??????????????????????not yet done */

/* -- */

    /* Parameter adjustments */
    ch -= 6;
    --h__;
    cv -= 4;
    --v;

    /* Function Body */
    ret_val = 1;
    w = h__[1];
    tl = h__[2];
    psi0 = h__[3];
    d0 = h__[4];
    z0 = h__[5];
/* -- track crosses r only, if */
/* -- d0<r and d0*w < 1 */
    if (abs(d0) < *r__) {
	srmd = sqrt(*r__ * *r__ - d0 * d0);
	sdtw = sqrt(1 - d0 * w);
	gamma = asin(w / 2. * srmd / sdtw) * 2.;
/* Computing 3rd power */
	d__1 = sdtw;
/* Computing 2nd power */
	d__2 = srmd * w / 2. / sdtw;
	dgamdu = .5 / srmd / (d__1 * (d__1 * d__1)) / sqrt(1. - d__2 * d__2);
	psi = psi0 + gamma;
	v[1] = xx;
	v[2] = yy;
	v[3] = zz;
	j[0] = 1.;
	j[3] = 0.;
	j[6] = 0.;
	j[9] = 0.;
	j[12] = 0.;
	j[1] = 0.;
	j[4] = 1.;
	j[7] = 0.;
	j[10] = 0.;
	j[13] = 0.;
/* Computing 3rd power */
	d__1 = sdtw;
/* Computing 2nd power */
	d__2 = srmd * w / 2. / sdtw;
	j[2] = (2. - d0 * w) * srmd / (d__1 * (d__1 * d__1) * 2. * sqrt(1. - 
		d__2 * d__2));
	j[5] = 0.;
	j[8] = 1.;
/* Computing 3rd power */
	d__1 = sdtw;
/* Computing 2nd power */
	d__2 = srmd * w / 2. / sdtw;
	j[11] = w * (d0 * -2. + d0 * d0 * w + *r__ * *r__ * w) / (srmd * 2. * 
		(d__1 * (d__1 * d__1)) * sqrt(1. - d__2 * d__2));
	j[14] = 1.;
	fvsabat_(&cv[4], j, &ch[6], &c__3, &c__5);
    } else {
	psi = psi0;
	for (i0 = 1; i0 <= 3; ++i0) {
	    for (i1 = 1; i1 <= 3; ++i1) {
		cv[i0 + i1 * 3] = ch[i0 + i1 * 5];
	    }
	}
    }
    return ret_val;
} /* fvrh2v_ */

integer fvhch_(doublereal *h__, doublereal *ch, doublereal *v, doublereal *q, 
	doublereal *c__, doublereal *d__, doublereal *e)
{
    /* System generated locals */
    integer ret_val;

    /* Local variables */
    extern integer fvh_(doublereal *, doublereal *, doublereal *), fvch_(
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static integer status;

/* -------------------------------------------------------------------------- */
/* -- calculate and return helix parameters h, Ch from v, q, C, D, E */

    /* Parameter adjustments */
    e -= 4;
    d__ -= 4;
    c__ -= 4;
    --q;
    --v;
    ch -= 6;
    --h__;

    /* Function Body */
    ret_val = 1;
    status = fvh_(&h__[1], &v[1], &q[1]);
    if ((status & 1) != 1) {
	ret_val = status;
    }
    status = fvch_(&ch[6], &v[1], &q[1], &c__[4], &d__[4], &e[4]);
    if ((status & 1) != 1) {
	ret_val = status;
    }
    return ret_val;
} /* fvhch_ */

integer fvh_(doublereal *h__, doublereal *v, doublereal *q)
{
    /* System generated locals */
    integer ret_val;

    /* Builtin functions */
    double sqrt(doublereal), atan2(doublereal, doublereal), cos(doublereal), 
	    sin(doublereal), atan(doublereal);

    /* Local variables */
    static doublereal r__, w, z__, tl, xi, xx, yy, phi, cxi, psi, sxi, oow, 
	    gamma;

/* -------------------------------------------------------------------------- */
/* -- calculate and return helix parameters h={w,tl,psi0,d0,z0} */
/* -- from v={vx,vy,vz} and q={w,tl,psi} */

    /* Parameter adjustments */
    --q;
    --v;
    --h__;

    /* Function Body */
    ret_val = 1;
    xx = v[1];
    yy = v[2];
    r__ = sqrt(xx * xx + yy * yy);
    phi = atan2(yy, xx);
    z__ = v[3];
    w = q[1];
    tl = q[2];
    psi = q[3];
    if (w != 0.) {
	oow = 1. / w;
	xi = psi - phi;
	cxi = cos(xi);
	sxi = sin(xi);
	gamma = atan(r__ * cxi / (oow - r__ * sxi));
	h__[1] = w;
	h__[2] = tl;
	h__[3] = psi - gamma;
	h__[4] = oow - (oow - r__ * sxi) / cos(gamma);
	h__[5] = z__ - gamma / w * tl;
    } else {
/* -- for neutral particles */
/* -- ?????? check this out */
	xi = psi - phi;
	sxi = sin(xi);
	h__[1] = w;
	h__[2] = tl;
	h__[3] = psi;
	h__[4] = r__ * sxi;
	h__[5] = z__;
    }
    return ret_val;
} /* fvh_ */

integer fvch_(doublereal *ch, doublereal *v, doublereal *q, doublereal *c__, 
	doublereal *d__, doublereal *e)
{
    /* System generated locals */
    integer ret_val;

    /* Local variables */
    static doublereal a[15]	/* was [5][3] */, b[15]	/* was [5][3] */, h0[
	    5];
    extern /* Subroutine */ int fvapb_(doublereal *, doublereal *, doublereal 
	    *, integer *, integer *);
    static doublereal temp55[25]	/* was [5][5] */;
    extern integer fvabh0_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    extern /* Subroutine */ int fvabct_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *), fvsabat_(
	    doublereal *, doublereal *, doublereal *, integer *, integer *), 
	    fvabtct_(doublereal *, doublereal *, doublereal *, doublereal *, 
	    integer *, integer *);

/* -------------------------------------------------------------------------- */
/* -- calculate and return cavariance matrix of helix parameters */
/* -- Ch = cov{Ch} from covariances of vertex v={vx,vy,vz} and momentum */
/* -- q={w,tl,psi}, */
/* -- C = cov{v}, D = cov{q}, E=cov{v,q} */

/* -- Ch = A.C.A^T + B.E.A^T + A.E^T.B^T + B.D.B^T */
    /* Parameter adjustments */
    e -= 4;
    d__ -= 4;
    c__ -= 4;
    --q;
    --v;
    ch -= 6;

    /* Function Body */
    ret_val = fvabh0_(a, b, h0, &v[1], &q[1]);
    fvsabat_(&ch[6], a, &c__[4], &c__5, &c__3);
    fvabct_(temp55, b, &e[4], a, &c__5, &c__3);
    fvapb_(&ch[6], &ch[6], temp55, &c__5, &c__5);
    fvabtct_(temp55, a, &e[4], b, &c__5, &c__3);
    fvapb_(&ch[6], &ch[6], temp55, &c__5, &c__5);
    fvsabat_(temp55, b, &d__[4], &c__5, &c__3);
    fvapb_(&ch[6], &ch[6], temp55, &c__5, &c__5);
    return ret_val;
} /* fvch_ */

integer fvabh0_(doublereal *a, doublereal *b, doublereal *h0, doublereal *v0, 
	doublereal *q0)
{
    /* System generated locals */
    integer ret_val;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal), atan2(doublereal, doublereal), d_mod(doublereal *
	    , doublereal *), cos(doublereal), sin(doublereal), atan(
	    doublereal);

    /* Local variables */
    static doublereal r__, w, z__, d0, z0, cg, sg, tl, xi, rw, xx, yy, phi, 
	    cxi, psi, sxi, oow, psi0, dgdw, dgdx, dgdy, drdx, drdy, gamma, 
	    dgdpsi, rdxidx, rdxidy, dgdvar0;

/* -------------------------------------------------------------------------- */
/* -- calculate and return A, B, h0 from v0, q0: */
/* -- */
/* -- calculate coefficients of taylor expansion for estimated */
/* -- track parameters p around point (v0,q0): */
/* -- */
/* --   p    =    h(v,q) + eps */
/* --     \approx h(v0,q0) + A.(v-v0) + B.(q-q0) + eps */
/* --        =    h0 + A.v + B.q + eps */
/* -- where */
/* --   A  =  dp/dv(v0,q0), */
/* --   B  =  dp/dq(v0,q0), and */
/* --   h0 =  h(v0,q0) - A.v0 - B.q0 are the */
/* -- derivatives of track parameters w/r to vertex position */
/* -- and track 3-momentum at space point v0 and estimated momentum q0 */
/* -- */

    /* Parameter adjustments */
    --q0;
    --v0;
    --h0;
    b -= 6;
    a -= 6;

    /* Function Body */
    ret_val = 1;
    xx = v0[1];
    yy = v0[2];
    r__ = sqrt(xx * xx + yy * yy);
    phi = atan2(yy, xx);
    z__ = v0[3];
    w = q0[1];
    tl = q0[2];
    psi = q0[3];
/* -- calculate some derived quantities */
    d__1 = psi - phi;
    xi = d_mod(&d__1, &c_b509);
    cxi = cos(xi);
    sxi = sin(xi);
    oow = 1 / w;
    rw = r__ * w;
    gamma = atan(r__ * cxi / (oow - r__ * sxi));
    sg = sin(gamma);
    cg = cos(gamma);
/* -- calc transformed quantities */
    psi0 = psi - gamma;
    d0 = oow - (oow - r__ * sxi) / cg;
    z0 = z__ - gamma / w * tl;
/* -- calc Jacobian */
    if (r__ != 0.) {
	drdx = xx / r__;
	drdy = yy / r__;
	rdxidx = yy / r__;
	rdxidy = -xx / r__;
    } else {
	drdx = 0.;
	drdy = 0.;
	rdxidx = 0.;
	rdxidy = 0.;
    }
    dgdvar0 = 1. / (rw * rw + 1. - rw * 2. * sxi);
    dgdx = dgdvar0 * (w * cxi * drdx + w * (rw - sxi) * rdxidx);
    dgdy = dgdvar0 * (w * cxi * drdy + w * (rw - sxi) * rdxidy);
    dgdw = dgdvar0 * r__ * cxi;
    dgdpsi = dgdvar0 * rw * (rw - sxi);
/* -- d w / d r, d phi, d z */
    a[6] = 0.;
    a[11] = 0.;
    a[16] = 0.;
/* -- d tl / d x, d y, d z */
    a[7] = 0.;
    a[12] = 0.;
    a[17] = 0.;
/* -- d psi0 / d x, d y, d z */
    a[8] = -dgdx;
    a[13] = -dgdy;
    a[18] = 0.;
/* -- d d0 / d x, d y, d z */
    a[9] = cxi / cg * rdxidx + sxi / cg * drdx - (oow - r__ * sxi) * sg / cg /
	     cg * dgdx;
    a[14] = cxi / cg * rdxidy + sxi / cg * drdy - (oow - r__ * sxi) * sg / cg 
	    / cg * dgdy;
    a[19] = 0.;
/* -- d z0 / d x, d y, d z */
    a[10] = -tl / w * dgdx;
    a[15] = -tl / w * dgdy;
    a[20] = 1.;
/* -- d w / d w, d tl, d psi */
    b[6] = 1.;
    b[11] = 0.;
    b[16] = 0.;
/* -- d tl / d w, d tl, d psi */
    b[7] = 0.;
    b[12] = 1.;
    b[17] = 0.;
/* -- d psi0 / d w, d tl, d psi */
    b[8] = -dgdw;
    b[13] = 0.;
    b[18] = 1. - dgdpsi;
/* -- d d0 / d w, d tl, d psi */
    b[9] = -oow * oow * (1. - 1. / cg) - (oow - r__ * sxi) / cg / cg * sg * 
	    dgdw;
    b[14] = 0.;
    b[19] = r__ * cxi / cg - (oow - r__ * sxi) / cg / cg * sg * dgdpsi;
/* -- d z0 / d w, d tl, d psi */
    b[10] = -tl / w * (dgdw - gamma / w);
    b[15] = -gamma / w;
    b[20] = -tl / w * dgdpsi;
    h0[1] = 0.;
    h0[2] = 0.;
    h0[3] = psi0 - a[8] * v0[1] - a[13] * v0[2] - b[8] * q0[1] - b[18] * q0[3]
	    ;
    h0[4] = d0 - a[9] * v0[1] - a[14] * v0[2] - b[9] * q0[1] - b[19] * q0[3];
    h0[5] = z0 - a[10] * v0[1] - a[15] * v0[2] - a[20] * v0[3] - b[10] * q0[1]
	     - b[15] * q0[2] - b[20] * q0[3];
    return ret_val;
} /* fvabh0_ */

integer fvcalcg_(doublereal *g, doublereal *cov, integer *n)
{
    /* System generated locals */
    integer g_dim1, g_offset, cov_dim1, cov_offset, ret_val;

    /* Local variables */
    extern integer fvsvdinv_(doublereal *, integer *);
    extern /* Subroutine */ int fvcopy_(doublereal *, doublereal *, integer *,
	     integer *);

/* ----------------------------------------------------------------------- */
/* -- calculate Gk = cov(pk)^(-1) for n x n matrix Cov */

/* -- invert covariance matrix */
    /* Parameter adjustments */
    cov_dim1 = *n;
    cov_offset = 1 + cov_dim1;
    cov -= cov_offset;
    g_dim1 = *n;
    g_offset = 1 + g_dim1;
    g -= g_offset;

    /* Function Body */
    fvcopy_(&g[g_offset], &cov[cov_offset], n, n);
    ret_val = fvsvdinv_(&g[g_offset], n);
    return ret_val;
} /* fvcalcg_ */

integer fvcalcw_(doublereal *w, doublereal *g, doublereal *b)
{
    /* System generated locals */
    integer ret_val;

    /* Local variables */
    extern integer fvsvdinv_(doublereal *, integer *);
    extern /* Subroutine */ int fvsatba_(doublereal *, doublereal *, 
	    doublereal *, integer *, integer *);

/* ----------------------------------------------------------------------- */
/* -- calculate W = (B^T.G.B)^(-1) */

    /* Parameter adjustments */
    b -= 6;
    g -= 6;
    w -= 4;

    /* Function Body */
    fvsatba_(&w[4], &b[6], &g[6], &c__5, &c__3);
    ret_val = fvsvdinv_(&w[4], &c__3);
    return ret_val;
} /* fvcalcw_ */

/* Subroutine */ int fvcalcgb_(doublereal *gb, doublereal *w, doublereal *b, 
	doublereal *g)
{
    extern /* Subroutine */ int fvamb_(doublereal *, doublereal *, doublereal 
	    *, integer *, integer *);
    static doublereal temp55[25]	/* was [5][5] */;
    extern /* Subroutine */ int fvsabat_(doublereal *, doublereal *, 
	    doublereal *, integer *, integer *);

/* ----------------------------------------------------------------------- */
/* -- calculate GB = G - G.B.W.B^T.G */

    /* Parameter adjustments */
    g -= 6;
    b -= 6;
    w -= 4;
    gb -= 6;

    /* Function Body */
    fvsabat_(temp55, &b[6], &w[4], &c__5, &c__3);
    fvsabat_(&gb[6], &g[6], temp55, &c__5, &c__5);
    fvamb_(&gb[6], &g[6], &gb[6], &c__5, &c__5);
    return 0;
} /* fvcalcgb_ */

/* Subroutine */ int fverror_(char *text, ftnlen text_len)
{
    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Fortran I/O blocks */
    static cilist io___287 = { 0, 6, 0, "(1x,a)", 0 };


/* ----------------------------------------------------------------------- */
/* -- output error string */

    s_wsfe(&io___287);
    do_fio(&c__1, text, text_len);
    e_wsfe();
    return 0;
} /* fverror_ */

/* Subroutine */ int fvsprt_(logical *p, integer *l)
{
/* ----------------------------------------------------------------------- */
/* -- set print flag and print unit */

    fvcom_1.print = *p;
    fvcom_1.plun = *l;
    return 0;
} /* fvsprt_ */

/* Subroutine */ int fvhelix2p4_(doublereal *t, doublereal *ct, doublereal *
	w2pt, doublereal *h__, doublereal *ch)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    double sin(doublereal), cos(doublereal);

    /* Local variables */
    static doublereal w, tl, ps, pt, px, py, pz, xy, cph, sph, sxx, sxy, syy, 
	    sxz, syz, szz, psi0, dpdk;

/* ----------------------------------------------------------------------- */
/* -- calculate 3-momentum from helix parameters parameters */
/* -- fill into 4-momentum array */
/* -- use w2pt for calculating pt from curvature w */

    /* Parameter adjustments */
    ch -= 6;
    --h__;
    ct -= 5;
    --t;

    /* Function Body */
    w = h__[1];
    tl = h__[2];
    psi0 = h__[3];
    sph = sin(psi0);
    cph = cos(psi0);
    pt = *w2pt / abs(w);
    px = pt * cph;
    py = pt * sph;
    pz = pt * tl;
/* -- calculate track momentum 3-vector and its error matrix */
    ps = *w2pt / w;
    dpdk = ps * ps / *w2pt;
    xy = ps * 2.f * dpdk * cph * sph * ch[16];
/* Computing 2nd power */
    d__1 = dpdk * cph;
/* Computing 2nd power */
    d__2 = ps * sph;
    sxx = d__1 * d__1 * ch[6] + d__2 * d__2 * ch[18] + xy;
    sxy = cph * sph * (dpdk * dpdk * ch[6] - ps * ps * ch[18]) + ps * dpdk * (
	    sph * sph - cph * cph) * ch[16];
/* Computing 2nd power */
    d__1 = dpdk * sph;
/* Computing 2nd power */
    d__2 = ps * cph;
    syy = d__1 * d__1 * ch[6] + d__2 * d__2 * ch[18] - xy;
    sxz = dpdk * dpdk * cph * tl * ch[6] - ps * dpdk * (cph * ch[11] - sph * 
	    tl * ch[16]) - ps * ps * sph * ch[17];
    syz = dpdk * dpdk * sph * tl * ch[6] - ps * dpdk * (sph * ch[11] + cph * 
	    tl * ch[16]) + ps * ps * cph * ch[17];
/* Computing 2nd power */
    d__1 = dpdk * tl;
    szz = d__1 * d__1 * ch[6] + ps * ps * ch[12] - ps * 2.f * dpdk * tl * ch[
	    11];
    t[1] = px;
    t[2] = py;
    t[3] = pz;
    ct[5] = sxx;
    ct[9] = sxy;
    ct[13] = sxz;
    ct[6] = sxy;
    ct[10] = syy;
    ct[14] = syz;
    ct[7] = sxz;
    ct[11] = syz;
    ct[15] = szz;
    return 0;
} /* fvhelix2p4_ */

/* Subroutine */ int fvhelix2p3_(doublereal *t, doublereal *ct, doublereal *
	w2pt, doublereal *h__, doublereal *ch)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    double sin(doublereal), cos(doublereal);

    /* Local variables */
    static doublereal w, tl, ps, pt, px, py, pz, xy, cph, sph, sxx, sxy, syy, 
	    sxz, syz, szz, psi0, dpdk;

/* ----------------------------------------------------------------------- */
/* -- calculate 3-momentum from helix parameters h, Ch */
/* -- fill into 3-momentum array */
/* -- use w2pt for calculating pt from curvature w */

    /* Parameter adjustments */
    ch -= 6;
    --h__;
    ct -= 4;
    --t;

    /* Function Body */
    w = h__[1];
    tl = h__[2];
    psi0 = h__[3];
    sph = sin(psi0);
    cph = cos(psi0);
    pt = *w2pt / abs(w);
    px = pt * cph;
    py = pt * sph;
    pz = pt * tl;
/* -- calculate track momentum 3-vector and its error matrix */
    ps = *w2pt / w;
    dpdk = ps * ps / *w2pt;
    xy = ps * 2.f * dpdk * cph * sph * ch[16];
/* Computing 2nd power */
    d__1 = dpdk * cph;
/* Computing 2nd power */
    d__2 = ps * sph;
    sxx = d__1 * d__1 * ch[6] + d__2 * d__2 * ch[18] + xy;
    sxy = cph * sph * (dpdk * dpdk * ch[6] - ps * ps * ch[18]) + ps * dpdk * (
	    sph * sph - cph * cph) * ch[16];
/* Computing 2nd power */
    d__1 = dpdk * sph;
/* Computing 2nd power */
    d__2 = ps * cph;
    syy = d__1 * d__1 * ch[6] + d__2 * d__2 * ch[18] - xy;
    sxz = dpdk * dpdk * cph * tl * ch[6] - ps * dpdk * (cph * ch[11] - sph * 
	    tl * ch[16]) - ps * ps * sph * ch[17];
    syz = dpdk * dpdk * sph * tl * ch[6] - ps * dpdk * (sph * ch[11] + cph * 
	    tl * ch[16]) + ps * ps * cph * ch[17];
/* Computing 2nd power */
    d__1 = dpdk * tl;
    szz = d__1 * d__1 * ch[6] + ps * ps * ch[12] - ps * 2.f * dpdk * tl * ch[
	    11];
    t[1] = px;
    t[2] = py;
    t[3] = pz;
    ct[4] = sxx;
    ct[7] = sxy;
    ct[10] = sxz;
    ct[5] = sxy;
    ct[8] = syy;
    ct[11] = syz;
    ct[6] = sxz;
    ct[9] = syz;
    ct[12] = szz;
    return 0;
} /* fvhelix2p3_ */

/* Subroutine */ int fvq2p4_(doublereal *p, doublereal *cp, doublereal *q, 
	doublereal *cq, doublereal *w2pt)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    double sin(doublereal), cos(doublereal);

    /* Local variables */
    static doublereal w, p0, tl, ps, pt, px, py, pz, xy, cp11, cp12, cp13, 
	    cp22, cp23, cp33, cph, sph, sxx, sxy, syy, sxz, syz, szz, dpdk;

/* -------------------------------------------------------------------------- */
/* -- calculate and return 3-momentum in four-vector p */
/* -- and its cov. matrix Cp in carth. coord. from */
/* -- q = {w, tl, psi} */
/* -- note: p, Cp is 4x4! */
/* -- */

    /* Parameter adjustments */
    cq -= 4;
    --q;
    cp -= 5;
    --p;

    /* Function Body */
    w = q[1];
    tl = q[2];
    p0 = q[3];
    sph = sin(p0);
    cph = cos(p0);
    pt = *w2pt / abs(w);
    cp11 = cq[4];
    cp12 = cq[7];
    cp13 = cq[10];
    cp22 = cq[8];
    cp23 = cq[11];
    cp33 = cq[12];
    px = pt * cph;
    py = pt * sph;
    pz = pt * tl;
/* -- calculate track momentum 3-vector and its error matrix */
    ps = *w2pt / w;
    dpdk = ps * ps / *w2pt;
    xy = ps * 2.f * dpdk * cph * sph * cp13;
/* Computing 2nd power */
    d__1 = dpdk * cph;
/* Computing 2nd power */
    d__2 = ps * sph;
    sxx = d__1 * d__1 * cp11 + d__2 * d__2 * cp33 + xy;
    sxy = cph * sph * (dpdk * dpdk * cp11 - ps * ps * cp33) + ps * dpdk * (
	    sph * sph - cph * cph) * cp13;
/* Computing 2nd power */
    d__1 = dpdk * sph;
/* Computing 2nd power */
    d__2 = ps * cph;
    syy = d__1 * d__1 * cp11 + d__2 * d__2 * cp33 - xy;
    sxz = dpdk * dpdk * cph * tl * cp11 - ps * dpdk * (cph * cp12 - sph * tl *
	     cp13) - ps * ps * sph * cp23;
    syz = dpdk * dpdk * sph * tl * cp11 - ps * dpdk * (sph * cp12 + cph * tl *
	     cp13) + ps * ps * cph * cp23;
/* Computing 2nd power */
    d__1 = dpdk * tl;
    szz = d__1 * d__1 * cp11 + ps * ps * cp22 - ps * 2.f * dpdk * tl * cp12;
    p[1] = px;
    p[2] = py;
    p[3] = pz;
    cp[5] = sxx;
    cp[9] = sxy;
    cp[13] = sxz;
    cp[6] = sxy;
    cp[10] = syy;
    cp[14] = syz;
    cp[7] = sxz;
    cp[11] = syz;
    cp[15] = szz;
    return 0;
} /* fvq2p4_ */

/* Subroutine */ int fvpe_(doublereal *p, doublereal *cp, integer *it, 
	doublereal *w2pt, doublereal *m, doublereal *ql, doublereal *cql)
{
    extern /* Subroutine */ int fvenergy_(doublereal *, doublereal *, 
	    doublereal *), fvq2p4_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);

/* ----------------------------------------------------------------------- */
/* -- calculate 4-momentum for track it */
    /* Parameter adjustments */
    cql -= 13;
    ql -= 4;
    cp -= 5;
    --p;

    /* Function Body */
    fvq2p4_(&p[1], &cp[5], &ql[*it * 3 + 1], &cql[(*it * 3 + 1) * 3 + 1], 
	    w2pt);
    fvenergy_(&p[1], &cp[5], m);
    return 0;
} /* fvpe_ */

/* Subroutine */ int fvenergy_(doublereal *t, doublereal *ct, doublereal *m)
{
    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal e, px, py, pz;

/* ----------------------------------------------------------------------- */
/* -- calculate energy from 3-momentum */
/* -- fill entries in t and Ct */

    /* Parameter adjustments */
    ct -= 5;
    --t;

    /* Function Body */
    px = t[1];
    py = t[2];
    pz = t[3];
    e = sqrt(*m * *m + px * px + py * py + pz * pz);
    t[4] = e;
    ct[17] = (px * ct[5] + py * ct[9] + pz * ct[13]) / e;
    ct[18] = (px * ct[9] + py * ct[10] + pz * ct[14]) / e;
    ct[19] = (px * ct[13] + py * ct[14] + pz * ct[15]) / e;
    ct[20] = (px * px * ct[5] + py * py * ct[10] + pz * pz * ct[15] + (px * (
	    py * ct[9] + pz * ct[13]) + py * pz * ct[14]) * 2.f) / e / e;
    ct[8] = ct[17];
    ct[12] = ct[18];
    ct[16] = ct[19];
    return 0;
} /* fvenergy_ */

/* Subroutine */ int fvq2p3_(doublereal *p, doublereal *cp, doublereal *q, 
	doublereal *cq, doublereal *w2pt)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    double sin(doublereal), cos(doublereal);

    /* Local variables */
    static doublereal w, p0, tl, ps, pt, px, py, pz, xy, cp11, cp12, cp13, 
	    cp22, cp23, cp33, cph, sph, sxx, sxy, syy, sxz, syz, szz, dpdk;

/* -------------------------------------------------------------------------- */
/* -- calculate and return 3-momentum vector p */
/* -- and its cov. matrix Cp in carth. coord. from */
/* -- q = {w, tl, psi} */
/* -- note: p, Cp is 3x3! */
/* -- */

    /* Parameter adjustments */
    cq -= 4;
    --q;
    cp -= 4;
    --p;

    /* Function Body */
    w = q[1];
    tl = q[2];
    p0 = q[3];
    sph = sin(p0);
    cph = cos(p0);
    pt = *w2pt / abs(w);
    cp11 = cq[4];
    cp12 = cq[7];
    cp13 = cq[10];
    cp22 = cq[8];
    cp23 = cq[11];
    cp33 = cq[12];
    px = pt * cph;
    py = pt * sph;
    pz = pt * tl;
/* -- calculate track momentum 3-vector and its error matrix */
    ps = *w2pt / w;
    dpdk = ps * ps / *w2pt;
    xy = ps * 2.f * dpdk * cph * sph * cp13;
/* Computing 2nd power */
    d__1 = dpdk * cph;
/* Computing 2nd power */
    d__2 = ps * sph;
    sxx = d__1 * d__1 * cp11 + d__2 * d__2 * cp33 + xy;
    sxy = cph * sph * (dpdk * dpdk * cp11 - ps * ps * cp33) + ps * dpdk * (
	    sph * sph - cph * cph) * cp13;
/* Computing 2nd power */
    d__1 = dpdk * sph;
/* Computing 2nd power */
    d__2 = ps * cph;
    syy = d__1 * d__1 * cp11 + d__2 * d__2 * cp33 - xy;
    sxz = dpdk * dpdk * cph * tl * cp11 - ps * dpdk * (cph * cp12 - sph * tl *
	     cp13) - ps * ps * sph * cp23;
    syz = dpdk * dpdk * sph * tl * cp11 - ps * dpdk * (sph * cp12 + cph * tl *
	     cp13) + ps * ps * cph * cp23;
/* Computing 2nd power */
    d__1 = dpdk * tl;
    szz = d__1 * d__1 * cp11 + ps * ps * cp22 - ps * 2.f * dpdk * tl * cp12;
    p[1] = px;
    p[2] = py;
    p[3] = pz;
    cp[4] = sxx;
    cp[7] = sxy;
    cp[10] = sxz;
    cp[5] = sxy;
    cp[8] = syy;
    cp[11] = syz;
    cp[6] = sxz;
    cp[9] = syz;
    cp[12] = szz;
    return 0;
} /* fvq2p3_ */

/* Subroutine */ int fvp3_(doublereal *p, doublereal *q, doublereal *w2pt)
{
    /* Builtin functions */
    double sin(doublereal), cos(doublereal);

    /* Local variables */
    static doublereal w, p0, tl, pt, cph, sph;

/* -------------------------------------------------------------------------- */
/* -- calculate and return 3-momentum vector p */
/* -- Cp in carth. coord. from */
/* -- q = {w, tl, psi} */
/* -- */

    /* Parameter adjustments */
    --q;
    --p;

    /* Function Body */
    w = q[1];
    tl = q[2];
    p0 = q[3];
    sph = sin(p0);
    cph = cos(p0);
    pt = *w2pt / abs(w);
    p[1] = pt * cph;
    p[2] = pt * sph;
    p[3] = pt * tl;
    return 0;
} /* fvp3_ */

/* Subroutine */ int fvq2pvec_(doublereal *p, doublereal *q, doublereal *w2pt)
{
    /* Builtin functions */
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static doublereal pt;

/* -------------------------------------------------------------------------- */
/* -- calculate and return 3-momentum vector p */
/* -- q = {w, tl, psi}, p = {px,py,pz} */
/* -- note: p, Cp is 3x3! */
/* -- */

    /* Parameter adjustments */
    --q;
    --p;

    /* Function Body */
    pt = *w2pt / abs(q[1]);
    p[1] = pt * cos(q[3]);
    p[2] = pt * sin(q[3]);
    p[3] = pt * q[2];
    return 0;
} /* fvq2pvec_ */

/* Subroutine */ int fvp2q_(doublereal *q, doublereal *p, doublereal *c__, 
	doublereal *w2pt)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), atan2(doublereal, doublereal);

    /* Local variables */
    static doublereal pt;

/* -------------------------------------------------------------------------- */
/* -- calculate q from 3-momentum vector p */
/* -- p = {px,py,pz}, q = {w, tl, psi} */
/* -- c is charge */
/* -- */

    /* Parameter adjustments */
    --p;
    --q;

    /* Function Body */
/* Computing 2nd power */
    d__1 = p[1];
/* Computing 2nd power */
    d__2 = p[2];
    pt = sqrt(d__1 * d__1 + d__2 * d__2);
/* -- ??????this is a kludge */
    if (*c__ != 0.) {
	q[1] = *c__ * *w2pt / pt;
    } else {
	q[1] = *w2pt / pt;
    }
    q[2] = p[3] / pt;
    q[3] = atan2(p[2], p[1]);
    return 0;
} /* fvp2q_ */

integer fvinv6s_(doublereal *c__, doublereal *d__, doublereal *e)
{
    /* System generated locals */
    integer ret_val;

    /* Local variables */
    extern integer fvsvdinv_(doublereal *, integer *);
    static doublereal m[36]	/* was [6][6] */;
    static integer i0, i1;

/* -------------------------------------------------------------------------- */
/* -- invert symmetric 6 x 6 matrix */
/* --       / C   E \ */
/* --   M = |  T     | */
/* --       \ E   D / */

/* -- fill up 6x6 matrix */
    /* Parameter adjustments */
    e -= 4;
    d__ -= 4;
    c__ -= 4;

    /* Function Body */
    for (i0 = 1; i0 <= 3; ++i0) {
	for (i1 = 1; i1 <= 3; ++i1) {
	    m[i0 + i1 * 6 - 7] = c__[i0 + i1 * 3];
	    m[i0 + 3 + (i1 + 3) * 6 - 7] = d__[i0 + i1 * 3];
	    m[i0 + (i1 + 3) * 6 - 7] = e[i0 + i1 * 3];
	    m[i1 + 3 + i0 * 6 - 7] = e[i0 + i1 * 3];
	}
    }
    ret_val = fvsvdinv_(m, &c__6);
/* -- get from 6x6 matric */
    for (i0 = 1; i0 <= 3; ++i0) {
	for (i1 = 1; i1 <= 3; ++i1) {
	    c__[i0 + i1 * 3] = m[i0 + i1 * 6 - 7];
	    d__[i0 + i1 * 3] = m[i0 + 3 + (i1 + 3) * 6 - 7];
	    e[i0 + i1 * 3] = m[i0 + (i1 + 3) * 6 - 7];
	}
    }
    return ret_val;
} /* fvinv6s_ */

/* C      integer function fvInv(A, n) */
/* CCCC-------------------------------------------------------------------------- */
/* CC -- invert symmetric n x n matrix A of order n, and return A^(-1) */
/* CC -- A will be overwritten! */
/* CC -- */
/* CC */
/* C      implicit none */
/* C      include 'fv.inc' */
/* C      integer n */
/* C      double precision A(n,n) */
/* C      integer ier */
/* C      integer temp(20) */
/* C */
/* C      fvInv = NORMAL */
/* CC#ifndef CRAY */
/* C      call DINV(n, A, n, temp, ier) */
/* CC      call DSINV(n, A, n, ier) */
/* CC#else */
/* CC      call RINV(n, A, n, temp, ier) */
/* CC#endif */
/* C      if (ier .ne. 0) then */
/* C         fvInv = ERROR */
/* C         return */
/* C      end if */
/* C      return */
/* C      end */
integer fvsvdinv_(doublereal *a, integer *n)
{

    /* System generated locals */
    integer a_dim1, a_offset, ret_val, i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);

    /* Local variables */
    static integer i__, j;
    static doublereal s, u[400]	/* was [20][20] */, v[400]	/* was [20][
	    20] */, w[20];
    static integer jj;
    static doublereal tmp[400]	/* was [20][20] */, wmax, jacob[20];
    extern integer fvsvd_(doublereal *, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *);
    static doublereal thresh;
    static integer status;

    /* Fortran I/O blocks */
    static cilist io___398 = { 0, 6, 0, 0, 0 };
    static cilist io___399 = { 0, 6, 0, 0, 0 };


/* -------------------------------------------------------------------------- */
/* -- invert symmetric n x n matrix A of order n, and return A^(-1) */
/* -- A will be overwritten! */
/* -- use a double precision version of the routine described in */
/* -- Press et al., Numerical Recipes, Cambridge Univ. Press */

    /* Parameter adjustments */
    a_dim1 = *n;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    ret_val = 1;
    if (*n > 20) {
	ret_val = 2;
	return ret_val;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	jacob[i__ - 1] = sqrt((d__1 = a[i__ + i__ * a_dim1], abs(d__1)));
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    u[i__ + j * 20 - 21] = a[i__ + j * a_dim1] / jacob[i__ - 1] / 
		    jacob[j - 1];
	}
    }
    status = fvsvd_(u, n, n, &c__20, &c__20, w, v);
    if (status != 1) {
	ret_val = status;
	return ret_val;
    }
/* -- edit singular values */
    wmax = 0.f;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (w[i__ - 1] > wmax) {
	    wmax = w[i__ - 1];
	}
    }
    thresh = wmax * 1e-14;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (w[i__ - 1] < thresh) {
	    s_wsle(&io___398);
	    do_lio(&c__9, &c__1, "fvSVDinv: editing singular values", (ftnlen)
		    33);
	    e_wsle();
	    s_wsle(&io___399);
	    do_lio(&c__9, &c__1, "----> ", (ftnlen)6);
	    do_lio(&c__5, &c__1, (char *)&w[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    e_wsle();
	    w[i__ - 1] = 0.;
	    ret_val = 3;
	}
    }
/* -- compose inverse matrix */
/* -- A <- A^(-1) = V.W^(-1).U^T */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    if (w[i__ - 1] != 0.) {
		tmp[i__ + j * 20 - 21] = u[j + i__ * 20 - 21] / w[i__ - 1];
	    } else {
		tmp[i__ + j * 20 - 21] = 0.;
	    }
	}
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    s = 0.;
	    i__3 = *n;
	    for (jj = 1; jj <= i__3; ++jj) {
		s += v[i__ + jj * 20 - 21] * tmp[jj + j * 20 - 21];
	    }
	    a[i__ + j * a_dim1] = s / jacob[i__ - 1] / jacob[j - 1];
	}
    }
    return ret_val;
} /* fvsvdinv_ */

integer fvluinv_(doublereal *y, doublereal *a, integer *n)
{
    /* System generated locals */
    integer y_dim1, y_offset, a_dim1, a_offset, ret_val, i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    extern /* Subroutine */ int fvlubksb_(doublereal *, integer *, integer *, 
	    doublereal *, doublereal *);
    extern integer fvludcmp_(doublereal *, integer *, integer *, doublereal *,
	     doublereal *);
    static doublereal d__;
    static integer i__, j;
    static doublereal tmp[36]	/* was [6][6] */, indx[6], jacob[6];
    static integer status;

/* -------------------------------------------------------------------------- */
/* -- invert n x n matrix A of order n, and return A^(-1) in Y */
/* -- using LU decomposition/backsubstitution */
/* -- A will be overwritten! */
/* -- use a double precision version of the routine described in */
/* -- Press et al., Numerical Recipes, Cambridge Univ. Press */

    /* Parameter adjustments */
    a_dim1 = *n;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    y_dim1 = *n;
    y_offset = 1 + y_dim1;
    y -= y_offset;

    /* Function Body */
    ret_val = 1;
    if (*n > 6) {
	ret_val = 2;
	return ret_val;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	jacob[i__ - 1] = sqrt((d__1 = a[i__ + i__ * a_dim1], abs(d__1)));
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    y[i__ + j * y_dim1] = 0.;
	    tmp[i__ + j * 6 - 7] = a[i__ + j * a_dim1] / jacob[i__ - 1] / 
		    jacob[j - 1];
	}
	y[i__ + i__ * y_dim1] = 1.;
    }
/* -- decompose the matrix just once */
    status = fvludcmp_(tmp, n, &c__6, indx, &d__);
    if ((status & 1) != 1) {
	ret_val = status;
	return ret_val;
    }
/* -- find invers by column */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	fvlubksb_(tmp, n, &c__6, indx, &y[j * y_dim1 + 1]);
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    y[i__ + j * y_dim1] = y[i__ + j * y_dim1] / jacob[i__ - 1] / 
		    jacob[j - 1];
	}
    }
    return ret_val;
} /* fvluinv_ */

/* Subroutine */ int fvcopy_(doublereal *a, doublereal *b, integer *n, 
	integer *m)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;

    /* Local variables */
    static integer i0, i1;

/* -------------------------------------------------------------------------- */
/* -- copy n x m matrix B into A */

/* #ifndef CRAY */
/*      call UCOPY(B, A, n*m*floatSize) */
/* #else */
    /* Parameter adjustments */
    b_dim1 = *n;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    a_dim1 = *n;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    i__1 = *n;
    for (i0 = 1; i0 <= i__1; ++i0) {
	i__2 = *m;
	for (i1 = 1; i1 <= i__2; ++i1) {
	    a[i0 + i1 * a_dim1] = b[i0 + i1 * b_dim1];
	}
    }
/* #endif */
    return 0;
} /* fvcopy_ */

/* Subroutine */ int fvzeroa_(doublereal *a, integer *n, integer *m)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j;

/* -------------------------------------------------------------------------- */
/* -- zero contents  into A for n x m matrices A */

    /* Parameter adjustments */
    a_dim1 = *n;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    a[i__ + j * a_dim1] = 0.;
	}
    }
    return 0;
} /* fvzeroa_ */

/* Subroutine */ int fvamb_(doublereal *x, doublereal *a, doublereal *b, 
	integer *n, integer *m)
{
    /* System generated locals */
    integer x_dim1, x_offset, a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j;

/* -------------------------------------------------------------------------- */
/* -- calculate X = A - B for n x m matrices X, A, B */
/* -- it is allowed for result matrix X to be A or B */

    /* Parameter adjustments */
    b_dim1 = *n;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    a_dim1 = *n;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    x_dim1 = *n;
    x_offset = 1 + x_dim1;
    x -= x_offset;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    x[i__ + j * x_dim1] = a[i__ + j * a_dim1] - b[i__ + j * b_dim1];
	}
    }
    return 0;
} /* fvamb_ */

/* Subroutine */ int fvnega_(doublereal *x, doublereal *a, integer *n, 
	integer *m)
{
    /* System generated locals */
    integer x_dim1, x_offset, a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j;

/* -------------------------------------------------------------------------- */
/* -- calculate X = -A for n x m matrices X, A */
/* -- it is allowed for result matrix X to be A */

    /* Parameter adjustments */
    a_dim1 = *n;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    x_dim1 = *n;
    x_offset = 1 + x_dim1;
    x -= x_offset;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    x[i__ + j * x_dim1] = -a[i__ + j * a_dim1];
	}
    }
    return 0;
} /* fvnega_ */

integer fvfta_(doublereal *x, doublereal *f, doublereal *a, integer *n, 
	integer *m)
{
    /* System generated locals */
    integer x_dim1, x_offset, a_dim1, a_offset, ret_val, i__1, i__2;

    /* Local variables */
    static integer i__, j;

/* -------------------------------------------------------------------------- */
/* -- calculate X = F*A for n x m matrices X, A and double precision number F */
/* -- it is allowed for result matrix X to be A */

    /* Parameter adjustments */
    a_dim1 = *n;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    x_dim1 = *n;
    x_offset = 1 + x_dim1;
    x -= x_offset;

    /* Function Body */
    ret_val = 1;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    x[i__ + j * x_dim1] = *f * a[i__ + j * a_dim1];
	}
    }
    return ret_val;
} /* fvfta_ */

/* Subroutine */ int fvapb_(doublereal *x, doublereal *a, doublereal *b, 
	integer *n, integer *m)
{
    /* System generated locals */
    integer x_dim1, x_offset, a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j;

/* -------------------------------------------------------------------------- */
/* -- calculate X = A + B for n x m matrices X, A, B */

    /* Parameter adjustments */
    b_dim1 = *n;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    a_dim1 = *n;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    x_dim1 = *n;
    x_offset = 1 + x_dim1;
    x -= x_offset;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    x[i__ + j * x_dim1] = a[i__ + j * a_dim1] + b[i__ + j * b_dim1];
	}
    }
    return 0;
} /* fvapb_ */

/* Subroutine */ int fvapbt_(doublereal *x, doublereal *a, doublereal *b, 
	integer *n, integer *m)
{
    /* System generated locals */
    integer x_dim1, x_offset, a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j;

/* -------------------------------------------------------------------------- */
/* -- calculate X = A + B^T for n x m matrices X, A, B */
/* -- note: X can be equal to A, but not to B */

    /* Parameter adjustments */
    b_dim1 = *n;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    a_dim1 = *n;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    x_dim1 = *n;
    x_offset = 1 + x_dim1;
    x -= x_offset;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    x[i__ + j * x_dim1] = a[i__ + j * a_dim1] + b[j + i__ * b_dim1];
	}
    }
    return 0;
} /* fvapbt_ */

/* Subroutine */ int fvab_(doublereal *x, doublereal *a, doublereal *b, 
	integer *n0, integer *n1, integer *n2)
{
    /* System generated locals */
    integer x_dim1, x_offset, a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, 
	    i__3;

    /* Local variables */
    static doublereal d__;
    static integer i__, j, k;

/* -------------------------------------------------------------------------- */
/* -- calculate X = A.B for n0 x n1 matrix A and n1 x n2 matrix B */
/* -- */
    /* Parameter adjustments */
    a_dim1 = *n0;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *n1;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    x_dim1 = *n0;
    x_offset = 1 + x_dim1;
    x -= x_offset;

    /* Function Body */
    i__1 = *n0;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n2;
	for (k = 1; k <= i__2; ++k) {
	    d__ = 0.;
	    i__3 = *n1;
	    for (j = 1; j <= i__3; ++j) {
		d__ += a[i__ + j * a_dim1] * b[j + k * b_dim1];
	    }
	    x[i__ + k * x_dim1] = d__;
	}
    }
    return 0;
} /* fvab_ */

/* Subroutine */ int fvatb_(doublereal *x, doublereal *a, doublereal *b, 
	integer *n0, integer *n1, integer *n2)
{
    /* System generated locals */
    integer x_dim1, x_offset, a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, 
	    i__3;

    /* Local variables */
    static doublereal d__;
    static integer i__, k, l;

/* -------------------------------------------------------------------------- */
/* -- calculate n1 x n2 matrix X = A^T.B */
/* -- with n0 x n1 matrix A and n0 x n2 matrix B */
/* -- */

    /* Parameter adjustments */
    a_dim1 = *n0;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *n0;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    x_dim1 = *n1;
    x_offset = 1 + x_dim1;
    x -= x_offset;

    /* Function Body */
    i__1 = *n1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n2;
	for (l = 1; l <= i__2; ++l) {
	    d__ = 0.;
	    i__3 = *n0;
	    for (k = 1; k <= i__3; ++k) {
		d__ += a[k + i__ * a_dim1] * b[k + l * b_dim1];
	    }
	    x[i__ + l * x_dim1] = d__;
	}
    }
    return 0;
} /* fvatb_ */

/* Subroutine */ int fvabt_(doublereal *x, doublereal *a, doublereal *b, 
	integer *n0, integer *n1, integer *n2)
{
    /* System generated locals */
    integer x_dim1, x_offset, a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, 
	    i__3;

    /* Local variables */
    static doublereal d__;
    static integer i__, k, l;

/* -------------------------------------------------------------------------- */
/* -- calculate n0 x n2 matrix X = A.B^T */
/* -- with n0 x n1 matrix A and n2 x n1 matrix B */
/* -- */

    /* Parameter adjustments */
    a_dim1 = *n0;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *n2;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    x_dim1 = *n0;
    x_offset = 1 + x_dim1;
    x -= x_offset;

    /* Function Body */
    i__1 = *n1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n2;
	for (l = 1; l <= i__2; ++l) {
	    d__ = 0.;
	    i__3 = *n0;
	    for (k = 1; k <= i__3; ++k) {
		d__ += a[i__ + k * a_dim1] * b[l + k * b_dim1];
	    }
	    x[i__ + l * x_dim1] = d__;
	}
    }
    return 0;
} /* fvabt_ */

/* Subroutine */ int fvatbc_(doublereal *x, doublereal *a, doublereal *b, 
	doublereal *c__, integer *n, integer *m)
{
    /* System generated locals */
    integer x_dim1, x_offset, a_dim1, a_offset, b_dim1, b_offset, c_dim1, 
	    c_offset, i__1, i__2, i__3;

    /* Local variables */
    static doublereal d__;
    static integer i__, j, k, l;
    static doublereal temp[400]	/* was [20][20] */;

/* -------------------------------------------------------------------------- */
/* -- calculate m x m matrix X = A^T.B.C */
/* -- for n x n matrix B */
/* -- and n x m matrices A,C */
/* -- */

    /* Parameter adjustments */
    b_dim1 = *n;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    c_dim1 = *n;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    a_dim1 = *n;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    x_dim1 = *m;
    x_offset = 1 + x_dim1;
    x -= x_offset;

    /* Function Body */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    d__ = 0.;
	    i__3 = *n;
	    for (l = 1; l <= i__3; ++l) {
		d__ += b[k + l * b_dim1] * c__[l + j * c_dim1];
	    }
	    temp[k + j * 20 - 21] = d__;
	}
    }
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    d__ = 0.;
	    i__3 = *n;
	    for (k = 1; k <= i__3; ++k) {
		d__ += a[k + i__ * a_dim1] * temp[k + j * 20 - 21];
	    }
	    x[i__ + j * x_dim1] = d__;
	}
    }
    return 0;
} /* fvatbc_ */

/* Subroutine */ int fvsatba_(doublereal *x, doublereal *a, doublereal *b, 
	integer *n, integer *m)
{
    /* System generated locals */
    integer x_dim1, x_offset, a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, 
	    i__3;

    /* Local variables */
    static doublereal d__;
    static integer i__, j, k, l;
    static doublereal temp[400]	/* was [20][20] */;

/* -------------------------------------------------------------------------- */
/* -- calculate m x m matrix X = A^T.B.A */
/* -- for symmetric n x n matrix B */
/* -- and n x m matrix A */
/* -- */

    /* Parameter adjustments */
    b_dim1 = *n;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    a_dim1 = *n;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    x_dim1 = *m;
    x_offset = 1 + x_dim1;
    x -= x_offset;

    /* Function Body */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    d__ = 0.;
	    i__3 = *n;
	    for (l = 1; l <= i__3; ++l) {
		d__ += b[k + l * b_dim1] * a[l + j * a_dim1];
	    }
	    temp[k + j * 20 - 21] = d__;
	}
    }
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__;
	for (j = 1; j <= i__2; ++j) {
	    d__ = 0.;
	    i__3 = *n;
	    for (k = 1; k <= i__3; ++k) {
		d__ += a[k + i__ * a_dim1] * temp[k + j * 20 - 21];
	    }
	    x[i__ + j * x_dim1] = d__;
	    if (i__ != j) {
		x[j + i__ * x_dim1] = d__;
	    }
	}
    }
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__ - 1;
	for (j = 1; j <= i__2; ++j) {
	    x[j + i__ * x_dim1] = x[i__ + j * x_dim1];
	}
    }
    return 0;
} /* fvsatba_ */

/* Subroutine */ int fvabct_(doublereal *x, doublereal *a, doublereal *b, 
	doublereal *c__, integer *m, integer *n)
{
    /* System generated locals */
    integer x_dim1, x_offset, a_dim1, a_offset, b_dim1, b_offset, c_dim1, 
	    c_offset, i__1, i__2, i__3;

    /* Local variables */
    static doublereal d__;
    static integer i__, j, k, l;
    static doublereal temp[400]	/* was [20][20] */;

/* -------------------------------------------------------------------------- */
/* -- calculate m x m matrix X = A.B.C^T */
/* -- for n x n matrix B */
/* -- and m x n matrices A, C */
/* -- */

    /* Parameter adjustments */
    x_dim1 = *m;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    c_dim1 = *m;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    b_dim1 = *n;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    a_dim1 = *m;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *m;
	for (l = 1; l <= i__2; ++l) {
	    d__ = 0.;
	    i__3 = *n;
	    for (j = 1; j <= i__3; ++j) {
		d__ += b[i__ + j * b_dim1] * c__[l + j * c_dim1];
	    }
	    temp[i__ + l * 20 - 21] = d__;
	}
    }
    i__1 = *m;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *m;
	for (l = 1; l <= i__2; ++l) {
	    d__ = 0.;
	    i__3 = *n;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		d__ += a[k + i__ * a_dim1] * temp[i__ + l * 20 - 21];
	    }
	    x[k + l * x_dim1] = d__;
	}
    }
    return 0;
} /* fvabct_ */

/* Subroutine */ int fvsabat_(doublereal *x, doublereal *a, doublereal *b, 
	integer *m, integer *n)
{
    /* System generated locals */
    integer x_dim1, x_offset, a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, 
	    i__3;

    /* Local variables */
    static doublereal d__;
    static integer i__, j, k, l;
    static doublereal temp[400]	/* was [20][20] */;

/* -------------------------------------------------------------------------- */
/* -- calculate m x m matrix X = A.B.A^T */
/* -- for a symmetric n x n matrix B and m x n matrix A */
/* -- */

    /* Parameter adjustments */
    x_dim1 = *m;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    b_dim1 = *n;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    a_dim1 = *m;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *m;
	for (l = 1; l <= i__2; ++l) {
	    d__ = 0.;
	    i__3 = *n;
	    for (j = 1; j <= i__3; ++j) {
		d__ += b[i__ + j * b_dim1] * a[l + j * a_dim1];
	    }
	    temp[i__ + l * 20 - 21] = d__;
	}
    }
    i__1 = *m;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k;
	for (l = 1; l <= i__2; ++l) {
	    d__ = 0.;
	    i__3 = *n;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		d__ += a[k + i__ * a_dim1] * temp[i__ + l * 20 - 21];
	    }
	    x[k + l * x_dim1] = d__;
	    if (l != k) {
		x[l + k * x_dim1] = d__;
	    }
	}
    }
    return 0;
} /* fvsabat_ */

/* Subroutine */ int fvabtct_(doublereal *x, doublereal *a, doublereal *b, 
	doublereal *c__, integer *m, integer *n)
{
    /* System generated locals */
    integer x_dim1, x_offset, a_dim1, a_offset, b_dim1, b_offset, c_dim1, 
	    c_offset, i__1, i__2, i__3;

    /* Local variables */
    static doublereal d__;
    static integer i__, j, k, l;
    static doublereal temp[400]	/* was [20][20] */;

/* -------------------------------------------------------------------------- */
/* -- calculate m x m matrix X = A.B^T.C^T */
/* -- for n x n matrix B */
/* -- and m x n matrices A, C */
/* -- */

    /* Parameter adjustments */
    x_dim1 = *m;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    c_dim1 = *m;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    b_dim1 = *n;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    a_dim1 = *m;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *m;
	for (l = 1; l <= i__2; ++l) {
	    d__ = 0.;
	    i__3 = *n;
	    for (j = 1; j <= i__3; ++j) {
		d__ += b[j + i__ * b_dim1] * c__[l + j * c_dim1];
	    }
	    temp[i__ + l * 20 - 21] = d__;
	}
    }
    i__1 = *m;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *m;
	for (l = 1; l <= i__2; ++l) {
	    d__ = 0.;
	    i__3 = *n;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		d__ += a[k + i__ * a_dim1] * temp[i__ + l * 20 - 21];
	    }
	    x[k + l * x_dim1] = d__;
	}
    }
    return 0;
} /* fvabtct_ */

doublereal fvprob_(doublereal *chi2, integer *ndof)
{
    /* System generated locals */
    real r__1;
    doublereal ret_val;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    extern doublereal prob_(real *, integer *);
    extern /* Subroutine */ int fverror_(char *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___466 = { 0, fvcom_1.text, 0, "(a,g12.4,a)", 132, 1 };


/* ---------------------------------------------------------------------- */
/* -- check for negative chi2 and p */

    if (*chi2 <= 0.) {
	s_wsfi(&io___466);
	do_fio(&c__1, "Error in fvProb: chi**2 = ", (ftnlen)26);
	do_fio(&c__1, (char *)&(*chi2), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, " <= 0.", (ftnlen)6);
	e_wsfi();
	fverror_(fvcom_1.text, (ftnlen)132);
	ret_val = 0.;
    } else {
	r__1 = (real) (*chi2);
	ret_val = prob_(&r__1, ndof);
    }
    return ret_val;
} /* fvprob_ */

doublereal fvpyth_(doublereal *a, doublereal *b)
{
    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal at, bt, ct;

/* -------------------------------------------------------------------- */
/* -- calculate sqrt(a*a + b*b) w/o destructive underflow or overflow */

    at = abs(*a);
    bt = abs(*b);
    if (at > bt) {
	ct = bt / at;
/* -- protect against underlow */
	if (ct < 1e-18) {
	    ct = 0.;
	}
	ret_val = at * sqrt(ct * ct + 1.);
    } else {
	if (bt != 0.) {
	    ct = at / bt;
	    if (ct < 1e-18) {
		ct = 0.;
	    }
	    ret_val = bt * sqrt(ct * ct + 1.);
	} else {
	    ret_val = 0.;
	}
    }
    return ret_val;
} /* fvpyth_ */

integer fvsvd_(doublereal *a, integer *m, integer *n, integer *mp, integer *
	np, doublereal *w, doublereal *v)
{
    /* System generated locals */
    integer a_dim1, a_offset, v_dim1, v_offset, ret_val, i__1, i__2, i__3;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double sqrt(doublereal), d_sign(doublereal *, doublereal *);

    /* Local variables */
    static doublereal c__, f, g, h__;
    static integer i__, j, k, l;
    static doublereal s, x, y, z__;
    static integer jj, nm;
    static doublereal rv1[100];
    static integer its;
    static doublereal scale, anorm;
    extern doublereal fvpyth_(doublereal *, doublereal *);

/* -------------------------------------------------------------------- */
/* -- do a singular value decomposition of matrix A */
/* -- */
/* -- origin: W.H. Press et al., Numerical Recipes, Cambridge Univ. Press */
/* -- */
/* -- Given a matrix A, with logical dimensions m by n */
/* -- and physical dimensions mp by np, this */
/* -- routine computes its singular value decomposition */
/* -- A=U.W.V^T. The matrix U replaces */
/* -- A on output. The diagonal matrix of singular values W */
/* -- is output as a vector W. The matrix */
/* -- V (not the transpose V^T) is output as V. */
/* -- m must be greater or equal to n; if it is smaller */
/* -- then A should be filled up to square with zero rows. */
/* -- */
/* -- maximum anticipated value of n */
    /* Parameter adjustments */
    v_dim1 = *np;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --w;
    a_dim1 = *mp;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    ret_val = 1;
    if (*m < *n) {
	ret_val = 2;
	return ret_val;
    }
/* -- Householder reduction to bidiagonal form */
    g = 0.;
    scale = 0.;
    anorm = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l = i__ + 1;
	rv1[i__ - 1] = scale * g;
	g = 0.;
	s = 0.;
	scale = 0.;
	if (i__ <= *m) {
	    i__2 = *m;
	    for (k = i__; k <= i__2; ++k) {
		scale += (d__1 = a[k + i__ * a_dim1], abs(d__1));
	    }
	    if (scale != 0.) {
		i__2 = *m;
		for (k = i__; k <= i__2; ++k) {
		    a[k + i__ * a_dim1] /= scale;
		    s += a[k + i__ * a_dim1] * a[k + i__ * a_dim1];
		}
		f = a[i__ + i__ * a_dim1];
		d__1 = sqrt(s);
		g = -d_sign(&d__1, &f);
		h__ = f * g - s;
		a[i__ + i__ * a_dim1] = f - g;
		if (i__ != *n) {
		    i__2 = *n;
		    for (j = l; j <= i__2; ++j) {
			s = 0.;
			i__3 = *m;
			for (k = i__; k <= i__3; ++k) {
			    s += a[k + i__ * a_dim1] * a[k + j * a_dim1];
			}
			f = s / h__;
			i__3 = *m;
			for (k = i__; k <= i__3; ++k) {
			    a[k + j * a_dim1] += f * a[k + i__ * a_dim1];
			}
		    }
		}
		i__2 = *m;
		for (k = i__; k <= i__2; ++k) {
		    a[k + i__ * a_dim1] = scale * a[k + i__ * a_dim1];
		}
	    }
	}
	w[i__] = scale * g;
	g = 0.;
	s = 0.;
	scale = 0.;
	if (i__ <= *m && i__ != *n) {
	    i__2 = *n;
	    for (k = l; k <= i__2; ++k) {
		scale += (d__1 = a[i__ + k * a_dim1], abs(d__1));
	    }
	    if (scale != 0.) {
		i__2 = *n;
		for (k = l; k <= i__2; ++k) {
		    a[i__ + k * a_dim1] /= scale;
		    s += a[i__ + k * a_dim1] * a[i__ + k * a_dim1];
		}
		f = a[i__ + l * a_dim1];
		d__1 = sqrt(s);
		g = -d_sign(&d__1, &f);
		h__ = f * g - s;
		a[i__ + l * a_dim1] = f - g;
		i__2 = *n;
		for (k = l; k <= i__2; ++k) {
		    rv1[k - 1] = a[i__ + k * a_dim1] / h__;
		}
		if (i__ != *m) {
		    i__2 = *m;
		    for (j = l; j <= i__2; ++j) {
			s = 0.;
			i__3 = *n;
			for (k = l; k <= i__3; ++k) {
			    s += a[j + k * a_dim1] * a[i__ + k * a_dim1];
			}
			i__3 = *n;
			for (k = l; k <= i__3; ++k) {
			    a[j + k * a_dim1] += s * rv1[k - 1];
			}
		    }
		}
		i__2 = *n;
		for (k = l; k <= i__2; ++k) {
		    a[i__ + k * a_dim1] = scale * a[i__ + k * a_dim1];
		}
	    }
	}
/* Computing MAX */
	d__3 = anorm, d__4 = (d__1 = w[i__], abs(d__1)) + (d__2 = rv1[i__ - 1]
		, abs(d__2));
	anorm = max(d__3,d__4);
    }
/* -- Accumulation of right-hand transformation */
    for (i__ = *n; i__ >= 1; --i__) {
	if (i__ < *n) {
	    if (g != 0.) {
		i__1 = *n;
		for (j = l; j <= i__1; ++j) {
/* -- double division to avoid possible underflow */
		    v[j + i__ * v_dim1] = a[i__ + j * a_dim1] / a[i__ + l * 
			    a_dim1] / g;
		}
		i__1 = *n;
		for (j = l; j <= i__1; ++j) {
		    s = 0.;
		    i__2 = *n;
		    for (k = l; k <= i__2; ++k) {
			s += a[i__ + k * a_dim1] * v[k + j * v_dim1];
		    }
		    i__2 = *n;
		    for (k = l; k <= i__2; ++k) {
			v[k + j * v_dim1] += s * v[k + i__ * v_dim1];
		    }
		}
	    }
	    i__1 = *n;
	    for (j = l; j <= i__1; ++j) {
		v[i__ + j * v_dim1] = 0.;
		v[j + i__ * v_dim1] = 0.;
	    }
	}
	v[i__ + i__ * v_dim1] = 1.;
	g = rv1[i__ - 1];
	l = i__;
    }
/* -- Accumulation of left-hand transformation */
    for (i__ = *n; i__ >= 1; --i__) {
	l = i__ + 1;
	g = w[i__];
	if (i__ < *n) {
	    i__1 = *n;
	    for (j = l; j <= i__1; ++j) {
		a[i__ + j * a_dim1] = 0.;
	    }
	}
	if (g != 0.) {
	    g = 1. / g;
	    if (i__ != *n) {
		i__1 = *n;
		for (j = l; j <= i__1; ++j) {
		    s = 0.;
		    i__2 = *m;
		    for (k = l; k <= i__2; ++k) {
			s += a[k + i__ * a_dim1] * a[k + j * a_dim1];
		    }
		    f = s / a[i__ + i__ * a_dim1] * g;
		    i__2 = *m;
		    for (k = i__; k <= i__2; ++k) {
			a[k + j * a_dim1] += f * a[k + i__ * a_dim1];
		    }
		}
	    }
	    i__1 = *m;
	    for (j = i__; j <= i__1; ++j) {
		a[j + i__ * a_dim1] *= g;
	    }
	} else {
	    i__1 = *m;
	    for (j = i__; j <= i__1; ++j) {
		a[j + i__ * a_dim1] = 0.;
	    }
	}
	a[i__ + i__ * a_dim1] += 1.;
    }
/* -- Diagonalization of the bidiagonal form */
    for (k = *n; k >= 1; --k) {
/* -- loop over singular values */
	for (its = 1; its <= 30; ++its) {
/* -- loop over allowed iterations */
	    for (l = k; l >= 1; --l) {
/* -- Test for splitting */
		nm = l - 1;
/* -- note tha rv1(1) is always zero */
		if ((d__1 = rv1[l - 1], abs(d__1)) + anorm == anorm) {
		    goto L2;
		}
		if ((d__1 = w[nm], abs(d__1)) + anorm == anorm) {
		    goto L1;
		}
	    }
L1:
	    c__ = 0.;
/* -- Cancellation of rv1(l) if l>1 */
	    s = 1.;
	    i__1 = k;
	    for (i__ = l; i__ <= i__1; ++i__) {
		f = s * rv1[i__ - 1];
		if (abs(f) + anorm != anorm) {
		    g = w[i__];
		    h__ = fvpyth_(&f, &g);
		    w[i__] = h__;
		    h__ = 1. / h__;
		    c__ = g * h__;
		    s = -(f * h__);
		    i__2 = *m;
		    for (j = 1; j <= i__2; ++j) {
			y = a[j + nm * a_dim1];
			z__ = a[j + i__ * a_dim1];
			a[j + nm * a_dim1] = y * c__ + z__ * s;
			a[j + i__ * a_dim1] = -(y * s) + z__ * c__;
		    }
		}
	    }
L2:
	    z__ = w[k];
	    if (l == k) {
/* -- Convergence */
		if (z__ < 0.) {
/* -- singular value is made nonnegative */
		    w[k] = -z__;
		    i__1 = *n;
		    for (j = 1; j <= i__1; ++j) {
			v[j + k * v_dim1] = -v[j + k * v_dim1];
		    }
		}
		goto L3;
	    }
	    if (its == 30) {
/* CCCC               call fvError('No convergence in 30 iterations') */
		ret_val = 2;
		return ret_val;
	    }
	    x = w[l];
/* -- shift from bottom 2-by-2 mirror */
	    nm = k - 1;
	    y = w[nm];
	    g = rv1[nm - 1];
	    h__ = rv1[k - 1];
	    f = ((y - z__) * (y + z__) + (g - h__) * (g + h__)) / (h__ * 2. * 
		    y);
	    g = fvpyth_(&f, &c_b134);
	    f = ((x - z__) * (x + z__) + h__ * (y / (f + d_sign(&g, &f)) - 
		    h__)) / x;
/* -- next QR transformation */
	    c__ = 1.;
	    s = 1.;
	    i__1 = nm;
	    for (j = l; j <= i__1; ++j) {
		i__ = j + 1;
		g = rv1[i__ - 1];
		y = w[i__];
		h__ = s * g;
		g = c__ * g;
		z__ = fvpyth_(&f, &h__);
		rv1[j - 1] = z__;
		c__ = f / z__;
		s = h__ / z__;
		f = x * c__ + g * s;
		g = -(x * s) + g * c__;
		h__ = y * s;
		y *= c__;
		i__2 = *n;
		for (jj = 1; jj <= i__2; ++jj) {
		    x = v[jj + j * v_dim1];
		    z__ = v[jj + i__ * v_dim1];
		    v[jj + j * v_dim1] = x * c__ + z__ * s;
		    v[jj + i__ * v_dim1] = -(x * s) + z__ * c__;
		}
		z__ = fvpyth_(&f, &h__);
		w[j] = z__;
/* -- rotation can be arbitrary if z=0 */
		if (z__ != 0.) {
		    z__ = 1. / z__;
		    c__ = f * z__;
		    s = h__ * z__;
		}
		f = c__ * g + s * y;
		x = -(s * g) + c__ * y;
		i__2 = *m;
		for (jj = 1; jj <= i__2; ++jj) {
		    y = a[jj + j * a_dim1];
		    z__ = a[jj + i__ * a_dim1];
		    a[jj + j * a_dim1] = y * c__ + z__ * s;
		    a[jj + i__ * a_dim1] = -(y * s) + z__ * c__;
		}
	    }
	    rv1[l - 1] = 0.;
	    rv1[k - 1] = f;
	    w[k] = x;
	}
L3:
	;
    }
    return ret_val;
} /* fvsvd_ */

integer fvludcmp_(doublereal *a, integer *n, integer *np, doublereal *indx, 
	doublereal *d__)
{
    /* System generated locals */
    integer a_dim1, a_offset, ret_val, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, k;
    static doublereal vv[100], dum, sum;
    static integer imax;
    static doublereal aamax;

/* ----------------------------------------------------------------------- */
/* -- do a LU-decomposition */
/* -- origin: W.H. Press et al., Numerical Recipes, Cambridge Univ. Press */
/* -- */
/* -- given an NxN matrix A, with physical dimensions NP, */
/* -- this routine replaces it by the LU */
/* -- decomposition of a rowwise permutation of itself. */
/* -- A and N are input. A is output, arranged */
/* -- as in equation (2.3.14) above; INDX is an output */
/* -- vector which records the row permutation */
/* -- effected by the partial pivoting; D is output */
/* -- as +- 1 depending on whatever the numer of */
/* -- row interchanges was even or odd, respectively. */
/* -- This routine is used in combination with */
/* -- LUBKSB to solve linear equations or invert a matrix. */

/* -- largest expected N, a small number */
/* -- VV stores the implicit scaling of each row */
    /* Parameter adjustments */
    --indx;
    a_dim1 = *np;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    ret_val = 1;
    *d__ = 1.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	aamax = 0.;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    if ((d__1 = a[i__ + j * a_dim1], abs(d__1)) > aamax) {
		aamax = (d__2 = a[i__ + j * a_dim1], abs(d__2));
	    }
	}
	if (aamax == 0.) {
/* -- matrix singular, no nonzero largest element */
	    ret_val = 2;
	    return ret_val;
	}
/* -- save the scaling */
	vv[i__ - 1] = 1. / aamax;
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* -- loop over columns of Crout's method */
	i__2 = j - 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* -- this is equation (2.3.12) except for i=j */
	    sum = a[i__ + j * a_dim1];
	    i__3 = i__ - 1;
	    for (k = 1; k <= i__3; ++k) {
		sum -= a[i__ + k * a_dim1] * a[k + j * a_dim1];
	    }
	    a[i__ + j * a_dim1] = sum;
	}
/* -- initialize for the search for the largest pivot element */
	aamax = 0.;
	i__2 = *n;
	for (i__ = j; i__ <= i__2; ++i__) {
/* -- this is i=j of equation (2.3.12) and i=j+1..n of equation (2.3.13) */
	    sum = a[i__ + j * a_dim1];
	    i__3 = j - 1;
	    for (k = 1; k <= i__3; ++k) {
		sum -= a[i__ + k * a_dim1] * a[k + j * a_dim1];
	    }
	    a[i__ + j * a_dim1] = sum;
/* -- figure of merit for the pivot */
	    dum = vv[i__ - 1] * abs(sum);
/* -- is it better than the best so far? */
	    if (dum >= aamax) {
		imax = i__;
		aamax = dum;
	    }
	}
/* -- do we need to interchange rows? */
	if (j != imax) {
/* -- yes, do so */
	    i__2 = *n;
	    for (k = 1; k <= i__2; ++k) {
		dum = a[imax + k * a_dim1];
		a[imax + k * a_dim1] = a[j + k * a_dim1];
		a[j + k * a_dim1] = dum;
	    }
/* -- change the parity of d */
	    *d__ = -(*d__);
/* -- also interchange the scale factor */
	    vv[imax - 1] = vv[j - 1];
	}
	indx[j] = (doublereal) imax;
/* -- now, finally, divide by the pivot element */
/* -- if the pivot element is zero the matrix is singulax */
/* -- (al leas to the precsision of the algorithm) */
/* -- for some applications on singular matrices it is desirable to */
/* -- substitute TINY for zero */
	if (a[j + j * a_dim1] == 0.) {
	    a[j + j * a_dim1] = 1e-30;
	}
	if (j != *n) {
	    dum = 1. / a[j + j * a_dim1];
	    i__2 = *n;
	    for (i__ = j + 1; i__ <= i__2; ++i__) {
		a[i__ + j * a_dim1] *= dum;
	    }
	}
/* -- go back to the next column in the reduction */
    }
    return ret_val;
} /* fvludcmp_ */

/* Subroutine */ int fvlubksb_(doublereal *a, integer *n, integer *np, 
	doublereal *indx, doublereal *b)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, ii, ll;
    static doublereal sum;

/* ----------------------------------------------------------------------- */
/* -- do a LU back substitution */
/* -- origin: W.H. Press et al., Numerical Recipes, Cambridge Univ. Press */
/* -- */
/* -- Solves the set of N linear equations A.X=B. */
/* -- Here A is input, not as the matrix A */
/* -- but rather as its LU decomposition, */
/* -- determined by the routine fvLUdcmp. INDX is input */
/* -- as the permutation vector returned by fvLUdcmp. */
/* -- B is input as the right-hand side vector */
/* -- B, and returns with the solution vector X. */
/* -- A, N, NP, and INDX are not modified by this */
/* -- routine and can be left in place for successive calls */
/* -- with different right-hand sides B. This */
/* -- routine takes into account the possibility that */
/* -- B will begin with many zero elements, so */
/* -- it is effiscient for use in matrix inversion. */

/* -- when ii is set to a positive value, it will become the index of the first */
/* -- nonvanishing element of b. We now do the forward substitution, */
/* -- equation (2.3.6). The only new wrinkle is to unscramble */
/* -- the permutation as we go */
    /* Parameter adjustments */
    --b;
    --indx;
    a_dim1 = *np;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    ii = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ll = (integer) indx[i__];
	sum = b[ll];
	b[ll] = b[i__];
	if (ii != 0) {
	    i__2 = i__ - 1;
	    for (j = ii; j <= i__2; ++j) {
		sum -= a[i__ + j * a_dim1] * b[j];
	    }
	} else if (sum != 0.) {
/* -- a nonzero element was encountered, so from now on we will have to */
/* -- do the sums in the loop above */
	    ii = i__;
	}
	b[i__] = sum;
    }
/* -- now we do the backsubstitution, equation 2.3.7 */
    for (i__ = *n; i__ >= 1; --i__) {
	sum = b[i__];
	if (i__ < *n) {
	    i__1 = *n;
	    for (j = i__ + 1; j <= i__1; ++j) {
		sum -= a[i__ + j * a_dim1] * b[j];
	    }
	}
/* -- store a component of the solution vector X */
	b[i__] = sum / a[i__ + i__ * a_dim1];
    }
/* -- all done! */
    return 0;
} /* fvlubksb_ */

/* Subroutine */ int fvtinv_(doublereal *a, doublereal *b, integer *n)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);

    /* Local variables */
    static doublereal c__[36]	/* was [6][6] */;
    static integer j, i0, i1;
    static doublereal sum;

    /* Fortran I/O blocks */
    static cilist io___506 = { 0, 6, 0, 0, 0 };
    static cilist io___507 = { 0, 6, 0, 0, 0 };
    static cilist io___508 = { 0, 6, 0, 0, 0 };


/* ------------------------------------------------------------------------- */
/* -- test inversion of a */
/* -- a = b^(-1) */
    /* Parameter adjustments */
    b_dim1 = *n;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    a_dim1 = *n;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    i__1 = *n;
    for (i0 = 1; i0 <= i__1; ++i0) {
	i__2 = *n;
	for (i1 = 1; i1 <= i__2; ++i1) {
	    sum = 0.;
	    i__3 = *n;
	    for (j = 1; j <= i__3; ++j) {
		sum += a[i0 + j * a_dim1] * b[j + i1 * b_dim1];
	    }
	    c__[i0 + i1 * 6 - 7] = sum;
	}
	if (a[i0 + i0 * a_dim1] < 0.) {
	    s_wsle(&io___506);
	    do_lio(&c__9, &c__1, "ERROR: error matrix a has neg. entries on "
		    "diag", (ftnlen)46);
	    e_wsle();
	}
	if (b[i0 + i0 * b_dim1] < 0.) {
	    s_wsle(&io___507);
	    do_lio(&c__9, &c__1, "ERROR: error matrix b has neg. entries on "
		    "diag", (ftnlen)46);
	    e_wsle();
	}
    }
    i__1 = *n;
    for (i0 = 1; i0 <= i__1; ++i0) {
	sum = 0.;
	i__2 = *n;
	for (i1 = 1; i1 <= i__2; ++i1) {
	    if (i0 != i1) {
		sum += a[i0 + i1 * a_dim1];
	    }
	}
	s_wsle(&io___508);
	do_lio(&c__3, &c__1, (char *)&i0, (ftnlen)sizeof(integer));
	do_lio(&c__5, &c__1, (char *)&a[i0 + i0 * a_dim1], (ftnlen)sizeof(
		doublereal));
	d__1 = sum / (doublereal) (*n);
	do_lio(&c__5, &c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	e_wsle();
    }
    return 0;
} /* fvtinv_ */

