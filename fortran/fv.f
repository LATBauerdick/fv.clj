C----------------------------------------------------------------------
C
C FV, package to fit vertices and find secondary vertices
C
C author: Lothar A.T. Bauerdick, CERN/PPE
C date:   Jan. 1991
C
C ref.:   R. Fruehwirt, ``Applications of Filter Methods...''
C         HEPHY-PUB 516/88, Vienna
C
C
C some notation:
C
C w:  omega, 1/R curvature >0 if helix turns anti-clockwise
C tl: tan(lambda), tangens of dip angle 
C psi: azimuth angle (in r-phi plane)
C psi0: psi @ point of closest approach to coordinate origin
C d0:   distance of helix at point of closest approach to origin in r-phi plane
C       sign convention of angular momentum Lz
C z0:   z coordinate of point of closest approach of helix to origin
C
C----------------------------------------------------------------------
      subroutine fv
      end
      block data FVBLKD
      include 'fv.inc'
      DATA 
     1	print/.false./,
     1	plun/41/
      end
      integer function fvFit(x,Cx,ql,Cql,chi2l,chi2t, 
     1                      nt,tList,x0,Cx0,hl,Ghl)
C-----------------------------------------------------------------------        
C -- calculate vertex x and its covariance matrix Cx (x = {x,y,z}),
C -- the list of momentum vectors ql and their covariance matrices
C -- Cql (q = {w, tl, psi}) and the total chi-square chi2t
C -- for nt tracks with numbers in tList
C -- from initial value for vertex x0, Cx0 and
C -- from list of track parameters hl and list of 
C -- inverse of covariance matrices of track parameters Ghl
C -- 3-momenta q are returned in {w, tl, psi}-system 
C --   (we do not know the transformation of w->pt at this point!)
C -- 
C -- input:
C --     nt             : number of tracks in tList
C --     tList({1..nt}) : index of tracks in hl, Ghl
C --     x0({1..3})     : initial guess of vertex position x0{x,y,z}
C --     Cx0({1..3},{1..3}) : covariance matrix of x0
C --     hl({1..5}, tList(i)) : track parameters 
C --                      p0={w,tl,psi0,d0,z0} of track i
C --     Ghl({1..5},{1..5}, tList(i)) : inverse of covariance matrix
C --                      (cov(p0(i)))^(-1) of p0 of track i
C -- output:
C --     x({1..3})              <- v({x,y,z}), coords of fitted vertx
C --     Cx({1..3})             <- cov{v} covariance matrix of v
C --     ql({1..3}, tList(i))   <- q({w,tl,psi})
C --     Cql({1..3},{1..3},tList(i)) <- covariance matrix of q
C --     chi2l(tList(i))        <- chi2 for this track belonging to v
C --     chi2t                  <- total chi2 for the vertex fit
C --     
C -- author: Lothar A.T. Bauerdick, CERN/PPE
C -- date:   Jan. 1991
C -- 
      implicit none
      include 'fv.inc'
      double precision x(dv),Cx(dv,dv),ql(dv,*),Cql(dv,dv,*),
     1                 chi2l(*),chi2t
      integer nt, tList(nt)
      double precision x0(dv), Cx0(dv,dv), hl(dh,*), Ghl(dh,dh,*)

      integer i0, i, it, status

      double precision v(dv), C(dv,dv), Gv(dv,dv),
     1                 v0(dv), C0(dv,dv), Gv0(dv,dv)
      double precision E(dv,dv), chi2
      double precision pp(dh),Cpp(dh,dh)
      double precision C5(dh,dh)

      integer fvq, fvhCh, fvLUinv
      external fvq, fvhCh, fvLUinv
      integer fvCalcG, fvFilter, fvSmooth, fvRemove
      external fvCalcG, fvFilter, fvSmooth, fvRemove
      double precision fvProb
      external fvProb

      fvFit = NORMAL
C
C -- the filter
C
C -- calculate start values v0, C0 from x0, Cx0
C -- v, Cv are in carth. coords too, ergo just copy
      call fvCopy(v0,x0,dv,1)
      call fvCopy(C0,Cx0,dv,dv)

      if (print) then
        write(plun,'(1x,a)') 
     1      'start values for Filter'
        write(plun,'(1x,3(g10.3,a,g10.3))') (v0(i0), ' +/-',
     1       sqrt(C0(i0,i0)) ,i0=1,dv)
      end if

C -- calculate inverse of covariance matrices for v0
        status = fvCalcG(Gv0,C0,dv)
        if (iand(status,1).NE.1) then
          fvFit = status
          return
        end if
        
      do i = 1, nt
        it = tList(i)
        if (it .NE. 0) then

          status = fvFilter(v,C,Gv,ql(1,it),Cql(1,1,it),E,chi2,
     1                      v0,Gv0,hl(1,it),Ghl(1,1,it))
          if (iand(status,1).NE.1) then
            write(text,'(a,i4)') 
     1        'Problems in fvFilter for Track ', it
            call fvError(text)
            fvFit = status
          end if

CLATB---------------------------
          if (print) then
             write(plun,'(1x,a,i4,a,g10.3)')
     1          'Filter result track: v and q for track', i,' chi2',chi2
             write(plun,'(1x,a,3(g10.3,a,g10.3))') "v0",(v0(i0), ' +/-',
     1       sqrt(C0(i0,i0)) ,i0=1,dv)
             status = fvLUinv(C5, Ghl(1,1,it),5)
            write(plun,'(1x,a,5(g10.3,a,g10.3))')"h ",(hl(i0,it),' +/-',
     1       sqrt(C5(i0,i0)) ,i0=1,dh)
             write(plun,'(1x,a,3(g10.3,a,g10.3))')"v ",(v(i0),' +/-',
     1       sqrt(C(i0,i0)) ,i0=1,dv)
            write(plun,'(1x,a,3(g10.3,a,g10.3))')"q ",(ql(i0,it),' +/-',
     1       sqrt(Cql(i0,i0,it)) ,i0=1,dv)
          end if
C -- use v and Gv as input vertex v0, Gv0 for next track
          call fvCopy(v0,v,dv,1)
          call fvCopy(Gv0,Gv,dv,dv)
          call fvCopy(C0,C,dv,dv)
        end if
      end do

CLATB print out filter end result

      write(plun,'(1x,a)')'-----------end Filter, start Smoother-------'


C
C -- the smoother
C -- v, C, Gv contain vertex, cov. matrix and cov^(-1)
      chi2t = 0d0
      do i = 1, nt
        it = tList(i)
        if (it .NE. 0) then
C -- some printout
          if (print) then
            write(plun,'(1x,a,i5)') 
     1        'Smoother for Track ', it
          end if

          status = fvSmooth(ql(1,it),Cql(1,1,it),E,chi2,
     1                      v,C,ql(1,it),hl(1,it),Ghl(1,1,it))
          if (iand(status,1).NE.1) then
            write(text,'(a,i4)') 
     1        'Problems in smoother for track ', it
            call fvError(text)
            fvFit = status
          end if

C -- calculate total chi2
          chi2t = chi2t + chi2

          if (print) then
C -- print out some results
            write(plun,'(1x,a,g10.3,a,g10.3)') 
     1        'chi2 =', chi2, 'prob (2 d.o.f.) =', fvProb(chi2, 2)
C -- get old track parameters error matrix
CCCCCCCCCC            call fvCopy(Cpp,Ghl(1,1,it),dh,dh)
            status = fvLUinv(Cpp, Ghl(1,1,it),dh)
            if (iand(status,1).NE.1) fvFit = status
            write(plun,'(1x,5(g10.3,a,g10.3))') (hl(i0,it), ' +/-',
     1         sqrt(Cpp(i0,i0)) ,i0=1,dh)
C -- calculate new track parameters and error matrix for printout
            status = fvhCh(pp,Cpp, v,ql(1,it),C,C ql(1,1,it),E)
            if (iand(status,1).NE.1) fvFit = status
            write(plun,'(1x,5(g10.3,a,g10.3))') (pp(i0), ' +/-',
     1       sqrt(Cpp(i0,i0)) ,i0=1,dh)
          end if
        end if
      end do

C -- save smoothed vertex position in x = {x, y, z}
      call fvCopy(x,v,dv,1)
      call fvCopy(Cx,C,dv,dv)

C -- smoothed vertex position still in v, C, Gv

      if (print) then
        write(plun,'(1x,a)') 
     1      'Smoothed result'
        write(plun,'(1x,a,g10.3,a,i3,a,i3,a,g10.3,a,g10.3,a)') 
     1      'chi2 =', chi2t, 
     1      'prob for ',2*nt,' (',2*nt-dv,') d.o.f. =',
     1      fvProb(chi2t,2*nt),
     1      ' (',fvProb(chi2t,2*nt-dv),')'
        write(plun,'(1x,3(g10.3,a,g10.3))') (v(i0), ' +/-',
     1     sqrt(C(i0,i0)) ,i0=1,dv)
      end if

C -- check for outliers
C -- calculate chi2 for each track to belong to this vertex
C -- set alpha = cut on probaility to zero, do not remove any track
      if (nt .le. 2) then
         chi2l(tList(1)) = 0d0
         chi2l(tList(2)) = 0d0
      else
        status = fvRemove(v0,C0,Gv0,ql,Cql,chi2l,nt,tList,
     1                    0d0,v,Gv,hl,Ghl)
        if (iand(status,1).NE.1) then
          fvFit = status
        end if
      end if

      return
      end
      integer function fvFitD(tau,sigTau,chi2, 
     1                        tau0,v0,C0,vp,Cp,p,Dp,Ep)
C-----------------------------------------------------------------------
C -- 
      implicit none
      include 'fv.inc'
      double precision tau,sigTau,chi2
      double precision tau0,v0(dv),C0(dv,dv),
     1                 vp(dv),Cp(dv,dv),p(dv),Dp(dv,dv),Ep(dv,dv)

      double precision pp, A(dv), v(dv), chi2Old, taup,
     1                 Jvp(dv,dv),Jp(dv,dv),C(dv,dv),G(dv,dv),
     1                 temp33(dv,dv),temp11(1)
      integer iter, status
      double precision chi2Cut
      parameter (chi2Cut = 0.1d0)
      integer maxIter
      parameter (maxIter=100)

      integer fvSVDinv
      external fvSVDinv


C -- fit tau for
C --    v = vp - tau*p/|p|
C -- where v is measured as v0 with covariance C0
C -- and vp,p are measured with Cp=cov{vp,vp},Dp=cov{p,p},Ep=cov{vp,p}

      fvFitD = NORMAL

      sigTau = 0d0

      pp = sqrt(p(1)**2+p(2)**2+p(3)**2)
C -- A = dv / dtau
      A(1) = -p(1)/pp
      A(2) = -p(2)/pp
      A(3) = -p(3)/pp

      tau = tau0
      iter = 0
      chi2Old = 1d10
  1   continue

C -- calculate jacobian J(v|vp,p) and transform covariance matrix
C -- for measurement of secondary vertex vp/p for given tau
      call fvJFitD(Jvp,Jp,p,pp,tau)

      call fvsABAT(C, Jvp,Cp,dv,dv)
      call fvsABAT(temp33, Jp,Dp,dv,dq)
      call fvAPB(C, C,temp33,dv,dv)

      call fvABCT(temp33, Jvp,Ep,Jp,dv,dv)
      call fvAPB(C, C,temp33,dv,dv)
      call fvAPBT(C, C,temp33,dv,dv)

C -- total covariance matrix is C + C0, invert it
      call fvAPB(G, C0,C,dv,dv)
      status = fvSVDinv(G,dv)
      if (status .eq. INFO) write(plun,'(1x,a)')
     1    'editing of singular values done during inv. of in tau fit'
      if (iand(status,1).NE.1) then
        fvFitD = status
        return
      end if

C -- calculate chi2

      v(1) = v0(1) - (vp(1) - tau*p(1)/pp)
      v(2) = v0(2) - (vp(2) - tau*p(2)/pp)
      v(3) = v0(3) - (vp(3) - tau*p(3)/pp)
      call fvsATBA(temp11, v,G,dv,1)
      chi2 = temp11(1)
      iter = iter+1

      if (abs(chi2Old - chi2) .LT. chi2Cut) then
C -- chi^2 does not change anymore, get out of here
        goto 99
      end if
      if (iter .GT. maxIter) then
C -- doesn't converge, return ERROR
        fvFitD = ERROR
        return
      end if
      chi2Old = chi2

C -- do the fit
C -- tau(i) = tau(i-1) + {A^T.G.A}^(-1).A^T.G.(v0-v(tau(i-1)))
C -- where A = dv/dtau | tau=tau(i-1)
C -- and G = (cov{v})^(-1)
      call fvsATBA(temp11, A,G,dv,1)
      taup = 1d0/temp11(1)

      call fvATBC(temp11, A,G,v,dv,1)

      taup = tau + taup*temp11(1)

      tau = taup
      goto 1

 99   continue

C -- re-calculate cov with most recent G
      call fvsATBA(temp11, A,G,dv,1)
      sigTau = sqrt(max(0d0,1d0/temp11(1)))

      return
      end
      subroutine fvJFitD(Jvp,Jp,p,pp,tau)
C-----------------------------------------------------------------------
C -- calculate Jacobian for v = vp - tau*p/|p|
C -- Jvp = J(v|vp), Jp = J(v|p)
C --
      implicit none
      include 'fv.inc'
      double precision Jvp(dv,dv),Jp(dv,dv),p(dv),pp,tau

      double precision pp2,tpp

      pp2 = pp**2
      tpp = tau/pp

      Jvp(1,1) = 1d0
      Jvp(1,2) = 0d0
      Jvp(1,3) = 0d0

      Jvp(2,1) = 0d0
      Jvp(2,2) = 1d0
      Jvp(2,3) = 0d0

      Jvp(3,1) = 0d0
      Jvp(3,2) = 0d0
      Jvp(3,3) = 1d0

      Jp(1,1) = -tpp*(1d0-p(1)*p(1)/pp2)
      Jp(1,2) = tpp*p(1)*p(2)/pp2
      Jp(1,3) = tpp*p(1)*p(3)/pp2

      Jp(2,1) = Jp(1,2)
      Jp(2,2) = -tpp*(1d0-p(2)*p(2)/pp2)
      Jp(2,3) = tpp*p(2)*p(3)/pp2

      Jp(3,1) = Jp(1,3)
      Jp(3,2) = Jp(2,3)
      Jp(3,3) = -tpp*(1d0-p(3)*p(3)/pp2)

      return
      end
      integer function fvDist(d,sigD, v0,Cv0,v1,Cv1)
C-----------------------------------------------------------------------
C -- calculate distance between two vertices
C
      implicit none
      include 'fv.inc'
      double precision d, sigD
      double precision v0(dv), Cv0(dv,dv), v1(dv), Cv1(dv)


      double precision dd(dv), temp11(1)

      fvDist = NORMAL

C -- distance
      d = sqrt((v0(1)-v1(1))**2 + 
     1         (v0(2)-v1(2))**2 + 
     1         (v0(3)-v1(3))**2)
C -- error
      dd(1) = (v0(1)-v1(1))/d
      dd(2) = (v0(2)-v1(2))/d
      dd(3) = (v0(3)-v1(3))/d
      call fvATBC(temp11,dd,Cv0,dd,3,1)
      sigD = temp11(1)
      call fvATBC(temp11,dd,Cv1,dd,3,1)
      sigD = sqrt(sigD + temp11(1))

      return
      end
      integer function fvInvM(m,sigm,
     1                        nt,tList,w2pt,ml,v,Cv,ql,Cql,Ghl)
C-----------------------------------------------------------------------
C -- return invariant mass of tracks
C
      implicit none
      include 'fv.inc'
      double precision m, sigm
      integer nt, tList(nt)
      double precision w2pt, ml(*), 
     1                 v(dv),Cv(dv,dv),
     1                 ql(dq,*),Cql(dq,dq,*),Ghl(dh,dh,*)


      double precision p(4), Cp(4,4), pi(4), Cpi(4,4)
      integer i,it, i0,i1

      fvInvM = NORMAL

      call fvZeroA(pi,4,1)
      call fvZeroA(Cpi,4,4)
      do i = 1, nt
        it = tList(i)
        call fvQ2P4(p(1),Cp(1,1), ql(1,it),Cql(1,1,it),w2pt)
        call fvEnergy(p(1), Cp(1,1), ml(it))
        do i0 = 1, 4
          pi(i0) = pi(i0) + p(i0)
          do i1 = 1, 4
            Cpi(i0,i1) = Cpi(i0,i1) + Cp(i0,i1)
          end do
        end do
      end do

C -- this does not take correlations between tracks into account
      m = pi(4)**2-pi(1)**2-pi(2)**2-pi(3)**2
      m = sqrt(max(m, 0d0))
      if (m .gt. 0d0) then
        sigm = 
     1      pi(1)*Cpi(1, 1)*pi(1) + 
     1	    pi(2)*Cpi(2, 2)*pi(2) + 
     1	    pi(3)*Cpi(3, 3)*pi(3) + 
     1	    pi(4)*Cpi(4, 4)*pi(4) +
     1	    2.*(pi(1)*(Cpi(1, 2)*pi(2) + 
     1            Cpi(1, 3)*pi(3) - 
     1            Cpi(1, 4)*pi(4)) +
     1          pi(2)*(Cpi(2, 3)*pi(3) - 
     1            Cpi(2, 4)*pi(4)) -
     1	        pi(3)*Cpi(3, 4)*pi(4))
        sigm = sqrt(max(sigm, 0d0))/m
      else
        sigm = 0d0
      end if


      return
      end
      integer function fvSumQ(p,Dp,Ep,
     1                        new,nt,tList,w2pt,v,C,ql,Cql,hl,Ghl)
C-----------------------------------------------------------------------
C -- sum up momenta of nt tracks in tList
C -- return summed 3-momentum and covariance Dp = cov{p} and Ep = cov{v,p}
C -- fill helix parameters and covariance matrix at position new
C -- in hl and Ghl
C --
C -- restrictions: assumes charge of +/- 1
C --               does not update ql,Cql
C --
      implicit none
      include 'fv.inc'      
      
      double precision p(dq),Dp(dq,dq),Ep(dq,dq)
      integer new,nt,tList(nt)
      double precision w2pt,v(dv),C(dv,dv),
     1                 ql(dq,*),Cql(dq,dq,*),hl(dh,*),Ghl(dh,dh,*)


      double precision 
     1                 D(dv,dv),E(dv,dv),
     1                 J0(dv,dv),J1(dv,dv),
     1                 q(dv),charge,Ch(dh,dh),
     1                 temp31(dv,1),temp33(dv,dv)

      integer i0,i1,it0,it1,status
      integer fvCovvq, fvCovqq, fvLUinv, fvhCh, fvJacobQ2P, fvh
      external fvCovvq, fvCovqq, fvLUinv, fvhCh, fvJacobQ2P, fvh

      fvSumQ = NORMAL

C -- check for overlap of tList(1..nt) and new
      do i0 = 1, nt
        if (new .EQ. tList(i0)) then
          fvSumQ = ERROR
        end if
      end do

C -- zero sum vector and error matrices
      call fvZeroA(p,dv,1)
      call fvZeroA(Dp,dv,dv)
      call fvZeroA(Ep,dv,dv)

C -- sum up charge, momentum and covariance matrix
      charge = 0d0
      do i0 = 1, nt
        it0 = tList(i0)
        charge = charge+sign(1d0,ql(1,it0))
        call fvQ2P3(temp31,temp33, ql(1,it0),Cql(1,1,it0),w2pt)
        call fvAPB(p,p,temp31,dv,1)
        call fvAPB(Dp,Dp,temp33,dv,dv)
      end do

C -- sum up covariance between different momenta
      do i0 = 1, nt
        it0 = tList(i0)
        status = fvJacobQ2P(J0,ql(1,it0),w2pt)
        if (iand(status,1).NE.1) fvSumQ = status
C -- covariance cov{v,q}
        status = fvCovvq(E, v,C,ql(1,it0),Ghl(1,1,it0))
        if (iand(status,1).NE.1) fvSumQ = status
        call fvABT(temp33, E,J0,dv,dv,dv)
        call fvAPB(Ep, Ep,temp33,dv,dv)

C -- get all Dij = cov{qi,qj}, i != j
        do i1 = i0+1,nt
          it1 = tList(i1)

C -- sum up correllation terms Dij = cov{qi,qj}:
C -- Dp = Dp + sum{i,j,I!=J}( Ji.Dij.Jj^T + Jj.Dij^T.Ji^T)
          status = fvCovqq(D, 
     1      v,C,ql(1,it0),Ghl(1,1,it0),ql(1,it1),Ghl(1,1,it1))
          if (iand(status,1).NE.1) fvSumQ = status
          status = fvJacobQ2P(J1,ql(1,it1),w2pt)
          if (iand(status,1).NE.1) fvSumQ = status
          call fvABCT(temp33, J0,D,J1,dv,dv)
          call fvAPB(Dp,Dp,temp33,dv,dv)
          call fvAPBT(Dp,Dp,temp33,dv,dv)

        end do
      end do

C -- should fill ql(1,new) and calculate
C -- D = dp/dq Dp dp/dq and fill into Cql(1,1,new)
C -- but we don't at the moment...
      call fvp2q(q, p,charge,w2pt)
      status = fvhCh(hl(1,new),Ch(1,1), v,q,C,Dp,Ep)
      if (iand(status,1).NE.1) fvSumQ = status

      status = fvLUinv(Ghl(1,1,new), Ch,dh)
      if (iand(status,1).NE.1) fvSumQ = status
CCCCC      call fvCopy(Ghl(1,1,new), Ch,dh,dh)

      return
      end
      integer function fvJacobQ2P(J,q,w2pt)
C-----------------------------------------------------------------------
C -- calculate Jacobian for q={w,tl,psi} -> p={px,py,pz} 
C -- i.e. J(P/Q)
C --
      implicit none
      include 'fv.inc'
      double precision J(dv,dv),q(dq),w2pt

      double precision w,tl,psi,cp,sp,pt,dptdw

      fvJacobQ2P = NORMAL

      w = q(1)
      tl = q(2)
      psi = q(3)
      cp = cos(psi)
      sp = sin(psi)

      pt = w2pt/abs(w)
      dptdw = -pt/w

      J(1,1) = dptdw*cp
      J(2,1) = dptdw*sp
      J(3,1) = dptdw*tl

      J(1,2) = 0d0
      J(2,2) = 0d0
      J(3,2) = pt

      J(1,3) = -pt*sp
      J(2,3) = pt*cp
      J(3,3) = 0d0

      return
      end
      integer function fvCovqq(Dij, v,C,qi,Gpi,qj,Gpj)
C-----------------------------------------------------------------------        
C -- calculate covariance matrix Dij between qi and qj
C --
      implicit none
      include 'fv.inc'
      double precision Dij(dq, dq)
      double precision v(dv), C(dv,dv),
     1                 qi(dq),Gpi(dh,dh),qj(dq),Gpj(dh,dh)


      double precision Ai(dh,dq),Bi(dh,dq),h0i(dh),
     1                 Aj(dh,dq),Bj(dh,dq),h0j(dh),
     1                 Wi(dq,dq),Wj(dq,dq),
     1                 temp55(dh,dh),tempp55(dh,dh),
     1                 temp33(dq,dq)

      integer status, fvABh0, fvCalcW
      external fvABh0, fvCalcW

      fvCovqq = NORMAL

C -- calculate coefficients for measurement equation, A, B, h0
      status = fvABh0(Ai,Bi,h0i, v,qi)
      if (iand(status,1).NE.1) fvCovqq = status
      status = fvABh0(Aj,Bj,h0j, v,qj)
      if (iand(status,1).NE.1) fvCovqq = status
C -- W
      status = fvCalcW(Wi,Gpi,Bi)
      if (status .eq. INFO) write(plun,'(1x,a)')
     1      'editing of singular values done during inv. of W'
      if (iand(status,1).NE.1) then
        fvCovqq = status
        return
      end if

      status = fvCalcW(Wj,Gpj,Bj)
      if (status .eq. INFO) write(plun,'(1x,a)')
     1      'editing of singular values done during inv. of W'
      if (iand(status,1).NE.1) then
        fvCovqq = status
        return
      end if

      call fvABCT(temp55, Ai,C,Aj,dh,dv)
      call fvABCT(tempp55, Gpi,temp55,Gpj,dh,dh)
      call fvATBC(temp33, Bi,tempp55,Bj,dh,dv)
      call fvABCT(Dij, Wi,temp33,Wj,dv,dv)

      return
      end
      integer function fvCovvq(E, v,C,q,Gh)
C-----------------------------------------------------------------------        
C -- calculate covariance matrix E between v and q 
C --
      implicit none
      include 'fv.inc'
      double precision E(dv,dv),
     1                 v(dv),C(dv,dv),q(dv),Gh(dh,dh)


      double precision A(dh,dv),B(dh,dv),h0(dh),
     1                 W(dv,dv), temp33(dv,dv)
      integer status
      integer fvABh0, fvCalcW
      external fvABh0, fvCalcW

      fvCovvq = NORMAL
C -- calculate coefficients for measurement equation, A, B, h0
      status = fvABh0(A,B,h0, v,q)
      if (iand(status,1).NE.1) fvCovvq = status
C -- W
      status = fvCalcW(W,Gh,B)
      if (status .eq. INFO) then
        write(plun,'(1x,a)')
     1      'editing of singular values done during inv. of W'
        fvCovvq = status
      end if
      if (iand(status,1).NE.1) then
        fvCovvq = status
        return
      end if

      call fvATBC(temp33, A,Gh,B,dh,dv)
      call fvABCT(E, C,temp33,W,dv,dv)
      call fvNegA(E,E,dv,dv)

      return
      end
      integer function fvFilter(v,C,Gv,q,D,E,chi2,
     1                          v0,Gv0,p0,Gp0)
C-----------------------------------------------------------------------        
C -- do the Kalman filter step for one track
C -- make use of fvFilterer
C
      implicit none
      include 'fv.inc'
      double precision v(dv),C(dv,dv),Gv(dv,dv),
     &  q(dq),D(dq,dv),E(dq,dq),chi2
      double precision v0(dv),Gv0(dv,dv),p0(dh),Gp0(dh,dh)


      double precision chi2Old,v00(dv),q00(dv),
     1                 A(dh,dv), B(dh,dv), h0(dh)
      integer status, iter

C -- iterate until chi2 change between iterations smaller chi2Cut
      double precision big, chi2Cut
      parameter (big = 1d10, chi2Cut = 0.5d0)

      integer fvq, fvABh0, fvFilterer, fvSVDfit
      external fvq, fvABh0, fvFilterer, fvSVDfit
C -- max #iteration in filter/smoother step
      integer maxFIter
      logical flag
      integer i0
      data maxFIter/10/
      data flag/.true./

      fvFilter = NORMAL

      call fvCopy(v00, v0,dv,1)
      status = fvq(q00, v00,p0)

      chi2 = big
      iter = 0
 1    continue
        chi2Old = chi2

C -- calculate coefficients for measurement equation, A, B, h0
        status = fvABh0(A,B,h0, v00,q00)
        if (iand(status,1).NE.1) fvFilter = status
        if (flag) then
           status = fvFilterer(v,C,Gv,q,D,E,chi2,
     1                         v0,Gv0,p0,Gp0,A,B,h0)
        else
           status = fvSVDfit(v,C,Gv,q,D,E,chi2, 
     1                       v0,Gv0,p0,Gp0,A,B,h0)
        end if
        if (iand(status,1).NE.1) fvFilter = status 

C -- use resultant vertex and 3-mom. for recalculating derivative
C -- and do it again
        call fvCopy(v00,v,dv,1)
        call fvCopy(q00,q,dv,1)

        iter = iter+1
C -- print out chi2 for each iteration
        if (print) then
          write(plun,'(1x,a,i5,a,g10.3)') 
     1      'Filter iteration ', iter, 
     1      ' yields chi2 ', chi2
        end if
      if (iter .LT. maxFIter .AND.
     1    abs(chi2Old-chi2) .GT. chi2Cut) goto 1


      return
      end
      integer function fvSmooth(q,D,E,chi2, v,C,q0,h,Gh)
C-----------------------------------------------------------------------        
C -- calculated smoothed momentum q and cov. matrix D,E 
C -- D = cov{q} and E = cov{v,q}
C -- with smoothed vertex position v,C
C -- use filtered momentum q0 = {w,tl,psi}
C -- return smoothed chi-square chi2
C
      implicit none
      include 'fv.inc'
      double precision q(dv),D(dq,dq),E(dv,dq),chi2
      double precision v(dv),C(dv,dv),q0(dq),h(dh),Gh(dh,dh)
      

      double precision A(dh,dv), B(dh,dv),h0(dh), chi2Old
      integer status

C -- iterate until chi2 change between iterations smaller chi2Cut
      double precision big, chi2Cut
      parameter (big = 1d10, chi2Cut = 0.5d0)
      integer iter
      integer fvABh0, fvSmoother
      external fvABh0, fvSmoother

C -- iterate only once at the moment
      integer maxSIter
      data maxSIter/1/

      fvSmooth = NORMAL

C -- use filtered q0 as expansion point for measurement equation
      call fvCopy(q, q0,dq,1)

      chi2 = big
      iter = 0
  1   continue
        chi2Old = chi2

C -- calculate coefficients for measurement equation, A, B, h0
        status = fvABh0(A,B,h0, v,q)
        if (iand(status,1).NE.1) fvSmooth = status
        status = fvSmoother(q,D,E,chi2,
     1                      v,C,h,Gh,A,B,h0)
        if (iand(status,1).NE.1) fvSmooth = status

C -- check on difference in chi2, ev. re-iterate 
C -- (maxSIter is 1 for the moment)
        iter = iter+1
      if (iter .LT. maxSIter .AND.
     1    abs(chi2Old-chi2) .GT. chi2Cut) goto 1

      return
      end
      integer function fvRemove(vp,Cvp,Gvp,ql,Cql,chi2l,nt,tList, 
     1                          alpha,v,Gv,hl,Ghl)
C-----------------------------------------------------------------------
C -- check for outliers and remove them from tList
C -- 
C -- calculate chi-square for each track to belong to vertex v
C -- remove tracks with chi-square probability larger than alpha
C -- return new vertex vp with covariance matrix Cvp, Gvp
C -- and updated 3-momenta ql, Cql with chi-square in chi2l
C -- removed tracks get a 0 in tList
C
      implicit none
      include 'fv.inc'
      double precision vp(dv),Cvp(dv,dv),Gvp(dv,dv),
     1                 ql(dq,*), Cql(dq,dq,*),chi2l(*)
      integer nt, tList(nt)
      double precision alpha
      double precision v(dv),Gv(dv,dv), 
     1                 hl(dh,*), Ghl(dh,dh,*)

      integer i, it, status, iMax
      double precision chi2Max

      double precision v0(dv),Gv0(dv,dv), 
     1                 A(dh,dv),B(dh,dv),h0(dh),chi2,E(dv,dv)

      logical removed

      integer fvABh0, fvRetlif, fvSmooth
      external fvABh0, fvRetlif, fvSmooth
      double precision prob
      double precision fvProb
      external fvProb

      fvRemove = NORMAL

      removed = .false.

C -- do not do remove with two tracks only
C -- set chi**2 to 1
      if (nt .le. 2) then
        it = tList(1)
        if (it .ne. 0) chi2l(it) = 2d0
        it = tList(2)
        if (it .ne. 0) chi2l(it) = 2d0
        call fvCopy(vp, v,dv,1)
        call fvCopy(Cvp,Gv,dv,dv)
        call fvCopy(Gvp,Gv,dv,dv)
        return
      end if

      call fvCopy(v0, v,dv,1)
      call fvCopy(Gv0,Gv,dv,dv)

  1   continue

      iMax = 0
      chi2Max = 0d0
      do i = 1, nt
        it = tList(i)
C -- check if this track has been removed
        if (it .NE. 0) then

C -- calculate coefficients for measurement equation, A, B, h0
C -- at smoothed vertex position v and smoothed 3-momentum q
          status = fvABh0(A,B,h0, v0,ql(1,it))
          if (iand(status,1).NE.1) then
            write(text,'(a,i4)') 
     1        'Problems in fvABh0 (outlier) for Track ', it
            call fvError(text)
            fvRemove = status
          end if
C -- calculate vertex vp, Cvp w/o track p0 and chi2
          status = fvRetlif(vp,Cvp,Gvp,chi2,
     1                 v0,Gv0,ql(1,it),hl(1,it),Ghl(1,1,it),A,B,h0)
          if (iand(status,1).NE.1) then
            write(text,'(a,i4)') 'Problems in fvRetlif for Track ', it
            call fvError(text)
            fvRemove = status
          end if
C -- save this chi2 in chi2l(it) and save i with max. chi2
          chi2l(it) = chi2
          if (chi2 .gt. chi2Max) then
            chi2Max = chi2
            iMax = i
          end if

          if (print) then
            prob = fvProb(chi2,2)
            write(plun,'(1x,a,i5,a,g10.3,a,g10.3)') 
     1        'Outlier Test for Track ', it, 
     1        ' chi2 ', chi2, ' prob. ', prob
          end if

        end if    !track has not been removed
      end do

C -- remove track with largest chi2 if prob(chi2) < alpha
      if (fvProb(chi2Max,2) .LT. alpha) then
        it = tList(iMax)
C -- remove track iMax from vertex and update all momenta
        status = fvRetlif(vp,Cvp,Gvp,chi2,
     1                  v0,Gv0,ql(1,it),hl(1,it),Ghl(1,1,it),A,B,h0)
        if (iand(status,1).NE.1) then
          write(text,'(a,i4)') 'Problems in fvRetlif for Track ', it
          call fvError(text)
          fvRemove = status
        end if
        tList(iMax) = 0
        do i = 1, nt
          it = tList(i)
          if (it .NE. 0) then
            status = fvSmooth(ql(1,it),Cql(1,1,it),E,chi2,
     1                        vp,Cvp,ql(1,it),hl(1,it),Ghl(1,1,it))
            if (iand(status,1).NE.1) then
              write(text,'(a,i4)') 
     1          'Problems in fvSmooth for Track ', it
              call fvError(text)
              fvRemove = status
            end if
          end if
        end do

C - do it again
        goto 1

      end if

      return
      end
      integer function fvSVDfit(v,C,Gv, q,D,E, chi2, 
     1                          v0,Gv0, h,Gh, A,B,h0)
C--------------------------------------------------------------------------
C -- calculate v, q, C, D, E, chi2 by solving the normal equation
C --
C -- / Gv0 + A^T.Gh.A   A^T.Gh.B \  / v \     / Gv0.v0 + A^T.Gh.(h-h0) \
C -- |                            | |    | =  |                         |
C -- \    B^T.Gh.A      B^T.Gh.B /  \ q /     \          B^T.Gh.(h-h0) /
C --
C -- Ansatz:
C --
C -- / v \    / C  E \    / Gv0.v0 + A^T.Gh.(h-h0) \
C -- |    | = |  T    | . |                         |
C -- \ q /    \ E  D /    \          B^T.Gh.(h-h0) /
C --
      implicit none
      include 'fv.inc'
      double precision v(dv),C(dv,dv),Gv(dv,dv),
     1                 q(dq),D(dq,dq),E(dv,dq),chi2,
     1                 v0(dv),Gv0(dv,dv), 
     1                 h(dh),Gh(dh,dh),
     1                 A(dh,dv),B(dh,dv),h0(dh)
      double precision 
     1                 m(dh), dm(dh), bv(dv), bq(dv),
     1                 temp33(dv,dv),
     1                 temp51(dh),
     1                 temp31(dv), tempp31(dv),
     1                 chi2ar(1,1)

      integer fvSVDinv, fvInv6s
      external fvSVDinv, fvInv6s
      integer status 

      fvSVDfit = NORMAL

C -- C <- Gv0 + A^T.Gh.A
      call fvsATBA(temp33, A,Gh,dh,dv)
      call fvAPB(C, Gv0,temp33,dv,dv)
C -- D <- B^T.Gh.B
      call fvsATBA(D, B,Gh,dh,dv)
C -- E <- A^T.Gh.B
      call fvATBC(E, A,Gh,B,dh,dv)

      status = fvInv6s(C,D,E)
      if (status .eq. INFO) write(plun,'(1x,a)')
     1      'editing of singular values done during inv. of CDE'
      if (iand(status,1).NE.1) then
        fvSVDfit = status
        return
      end if

C -- Gv = C^(-1)
      call fvCopy(Gv, C,dv,dv)
      status = fvSVDinv(Gv,dv)
      if (status .eq. INFO) write(plun,'(1x,a)')
     1    'editing of singular values done inverting C in fvSVDfit'
      if (iand(status,1).NE.1) then
        fvSVDfit = status
      end if

C -- m = h - h0
      call fvAMB(m, h,h0,dh,1)

C -- bv = Gv0.v0 + A^T.Gh.(h - h0)
C -- bq = B^T.Gh.(h - h0)
      call fvAB(temp31,Gv0,v0,dv,dv,1)
      call fvAB(temp51,Gh,m,dh,dh,1)
      call fvATB(tempp31,A,temp51,dh,dv,1)
      call fvAPB(bv,temp31,tempp31,dv,1)
      call fvATB(bq,B,temp51,dh,dv,1)

C -- v = C.bv + E.bq
      call fvAB(v,C,bv,dv,dv,1)
      call fvAB(temp31,E,bq,dv,dv,1)
      call fvAPB(v,v,temp31,dv,1)

C -- q = E^T.bv + D.bq
      call fvATB(q,E,bv,dv,dv,1)
      call fvAB(temp31,D,bq,dv,dv,1)
      call fvAPB(q,q,temp31,dv,1)

C -- dm = h - h0 - A.v
      call fvAB(temp51, A,v,dh,dv,1)
      call fvAMB(dm, m,temp51,dh,1)

C -- chi2 = (h - h0 - A.v - B.q)^T.Gh.(h - h0 - A.v - B.q) + 
C --        (v-v0)^T.Gv0.(v-v0)
      call fvAB(temp51, B,q,dh,dv,1)
      call fvAMB(temp51, dm,temp51,dh,1)
      call fvsATBA(chi2ar, temp51,Gh,dh,1)
      chi2 = chi2ar(1,1)

      call fvAMB(temp31, v,v0,dv,1)
      call fvsATBA(chi2ar, temp31,Gv0,dv,1)
      chi2 = chi2 + chi2ar(1,1)

      return
      end
      integer function fvFilterer(v,C,Gv,q,D,E,chi2, 
     1                            v0,Gv0,p0,G,A,B,h0)
C--------------------------------------------------------------------------
C -- Kalman filter for vertex fit
C -- 
C -- needs 2 inversions of 3x3-matrices
C --
C -- author: Lothar A.T. Bauerdick, CERN/PPE
C -- date:   Jan. 1991
C --
      implicit none
      include 'fv.inc'
      integer i0, i1
      double precision v(dv), C(dv,dv), Gv(dv,dv), 
     1                 q(dq), D(dq,dq), E(dv,dq), chi2
      double precision v0(dv), Gv0(dv,dv), 
     1                 p0(dh), h0(dh), G(dh,dh),
     1                 A(dh,dv), B(dh,dq)
      double precision W(dv,dv), GB(dh,dh), 
     1                 m(dh), dm(dh),
     1                 temp55(dh,dh), tempp55(dh,dh),
     1                 temp33(dv,dv),
     1                 temp51(dh),
     1                 temp31(dv), tempp31(dv),
     1                 chi2ar(1,1)

      integer fvSVDinv
      external fvSVDinv
      integer status 
      integer fvCalcW
      external fvCalcW

      fvFilterer = NORMAL

C -- W
      status = fvCalcW(W,G,B)
      if (status .eq. INFO) then
          write(text,'(a)')
     1      'editing of singular values done during inv. of W'
          call fvError(text)
C -- ????
C -- for now: forget about the fit...
          fvFilterer = ERROR
          return
      end if
      if (iand(status,1).NE.1) then
        fvFilterer = status
        return
      end if

C -- GB
      call fvCalcGB(GB, W,B,G)
C -- C
      call fvsATBA(Gv, A,GB,dh,dv)
      call fvAPB(Gv, Gv0,Gv,dv,dv)
      call fvCopy(C, Gv,dv,dv)
      status = fvSVDinv(C, dv)
      if (status .eq. INFO) write(plun,'(1x,a)')
     1  'editing of singular values done during inv. of C^(-1)'
      if (iand(status,1).NE.1) then
        fvFilterer = status
        return
      end if
C -- m = p0 - h0
      call fvAMB(m,p0,h0,dh,1)
C -- v
      call fvAB(temp51,GB,m,dh,dh,1)
      call fvATB(temp31,A,temp51,dh,dv,1)
      call fvAB(tempp31,Gv0,v0,dv,dv,1)
      call fvAPB(temp31,tempp31,temp31,dv,1)

      call fvAB(v,C,temp31,dv,dv,1)

C -- dm = p0 - h0 - A.v
      call fvAB(temp51,A,v,dh,dv,1)
      call fvAMB(dm,m,temp51,dh,1)
C -- q
      call fvAB(temp51,G,dm,dh,dh,1)
      call fvATB(temp31,B,temp51,dh,dv,1)
      call fvAB(q,W,temp31,dv,dv,1)

C -- D
      call fvsABAT(temp55,A,C,dh,dv)
      call fvsATBA(tempp55, G,temp55,dh,dh)
      call fvsATBA(temp33, B,tempp55,dh,dv)
      call fvsATBA(D, W,temp33,dv,dv)
      call fvAPB(D, W,D,dv,dv)

C -- E
      call fvATBC(temp33, B,G,A,dh,dv)
      call fvATBC(E,W,temp33,C,dv,dv)
      call fvNegA(E,E,dv,dv)

C -- chi2
      call fvAB(temp51, B,q,dh,dv,1)
      call fvAMB(temp51, dm, temp51, dh, 1)
      call fvsATBA(chi2ar, temp51,G,dh,1)
      chi2 = chi2ar(1,1)

      call fvAMB(temp31, v,v0,dv,1)
      call fvsATBA(chi2ar, temp31,Gv0,dv,1)
      chi2 = chi2 + chi2ar(1,1)
      
CLATB print A B h0 W GB dm C D E
C     write(plun,'(1x,a,3(g10.3))') "v0",v0
C     write(plun,'(1x,a,3(g10.3))') "Gv0",Gv0
C     write(plun,'(1x,a,5(g10.3))') "A ",A
C     write(plun,'(1x,a,5(g10.3))') "B ",B
C     write(plun,'(1x,a,5(g10.3))') "h0 ",h0
C     write(plun,'(1x,a,3(g10.3))') "W ",W
C     write(plun,'(1x,a,5(g10.3))') "GB ",GB
C     write(plun,'(1x,a,5(g10.3))') "dm ",dm
C     write(plun,'(1x,a,3(g10.3))') "D ",D
C     write(plun,'(1x,a,3(g10.3))') "E ",E
C     write(plun,'(1x,a,3(g10.3))') "v ",v
C     write(plun,'(1x,a,3(g10.3))') "C ",C
C     write(plun,'(1x,a,3(g10.3))') "chi2",chi2
      return
      end
      integer function fvSmoother(q, D, E, chi2, 
     1                            v,C,h,G,A,B,h0)
C--------------------------------------------------------------------------
C -- smoother
C --
C -- author: Lothar A.T. Bauerdick, CERN/PPE
C -- date:   Jan. 1991
C --
      implicit none
      include 'fv.inc'
      double precision q(dv), D(dv,dv), E(dv,dq), chi2,
     1                 v(dv), C(dv,dv), h(dh), G(dh,dh),
     1                 A(dh,dv), B(dh,dq), h0(dh)
      double precision W(dv,dv), GB(dh,dh), 
     1                 m(dh), dm(dh),
     1                 temp55(dh,dh), tempp55(dh,dh),
     1                 temp33(dv,dv),
     1                 temp51(dh),
     1                 temp31(dv),
     1                 chi2ar(1,1),
     1                 Gh(dh,dh)

      integer fvSVDinv
      external fvSVDinv
      integer status 
      integer fvCalcG, fvCalcW, fvCh, fvh
      external fvCalcG, fvCalcW, fvCh, fvh

      fvSmoother = NORMAL

C -- W
      status = fvCalcW(W,G,B)
      if (status .eq. INFO) write(plun,'(1x,a)')
     1      'editing of singular values done during inv. of W'
      if (iand(status,1).NE.1) then
        fvSmoother = status
        return
      end if
C -- GB
      call fvCalcGB(GB,W,B,G)
C -- m = h - h0
      call fvAMB(m, h,h0,dh,1)
C -- dm = h - h0 - A.v
      call fvAB(temp51,A,v,dh,dv,1)
      call fvAMB(dm,m,temp51,dh,1)
C -- q
      call fvAB(temp51,G,dm,dh,dh,1)
      call fvATB(temp31,B,temp51,dh,dv,1)
      call fvAB(q,W,temp31,dv,dv,1)
C -- D
      call fvsABAT(temp55, A,C,dh,dv)
      call fvsATBA(tempp55, G,temp55,dh,dh)
      call fvsATBA(temp33, B,tempp55,dh,dv)
      call fvsATBA(D, W,temp33,dv,dv)
      call fvAPB(D,W,D,dv,dv)

C -- E
      call fvATBC(temp33,B,G,A,dh,dv)
      call fvATBC(E,W,temp33,C,dv,dv)
      call fvNegA(E,E,dv,dv)

C -- chi2 = (h - h(v,q))^T.Gh(v,q).(h - h(v,q))
C -- where we calculate Gh(v,q) from the fit-result covariance 
C -- matrices C,D,E for the smoothed results v,q
CCCCC      call fvAB(temp51, B,q,dh,dv,1)
CCCCC      call fvAMB(temp51, dm,temp51,dh,1)

CCC      status = fvCh(Gh, v,q,C,D,E)
CCC      if (status .NE. NORMAL) then
CCC         fvSmoother = status
CCC      endif
CCC      status = fvSVDinv(Gh,dh)
CCC      if (status .NE. NORMAL) then
CCC         fvSmoother = status
CCC      endif

      status = fvh(temp51, v,q)
      if (iand(status,1).NE.1) fvSmoother = status
      call fvAMB(temp51, h,temp51,dh,1)
CCC      call fvsATBA(chi2ar, temp51,Gh,dh,1)
      call fvsATBA(chi2ar, temp51,G,dh,1)
      chi2 = chi2ar(1,1)

CLATB print A B h0 W GB dm C D E
      write(plun,'(1x,a,3(g10.3))') "v",v
      write(plun,'(1x,a,3(g10.3))') "C",C
      write(plun,'(1x,a,5(g10.3))') "A ",A
      write(plun,'(1x,a,5(g10.3))') "B ",B
      write(plun,'(1x,a,5(g10.3))') "h0 ",h0
      write(plun,'(1x,a,3(g10.3))') "W ",W
      write(plun,'(1x,a,5(g10.3))') "dm ",dm
      write(plun,'(1x,a,3(g10.3))') "E ",E
      write(plun,'(1x,a,3(g10.3))') "q ",q
      write(plun,'(1x,a,3(g10.3))') "D ",D
      write(plun,'(1x,a,3(g10.3))') "chi2",chi2
      return
      end
      integer function fvRetlif(vp, Cp, Gvp, chi2, 
     1                          v, Gv, q, p0, G, A, B, h0)
C--------------------------------------------------------------------------
C -- inverse Kalman filter, removes track p0, G from vertex v, Gv
C -- returns new vertex vp, Cp and the chi2 for track belonging to
C -- new vertex
C --
C -- author: Lothar A.T. Bauerdick, CERN/PPE
C -- date:   Jan. 1991
C --
      implicit none
      include 'fv.inc'
      double precision vp(dv), Cp(dv,dv), Gvp(dv,dv), chi2,
     1                 v(dv), Gv(dv,dv), 
     1                 q(dq), p0(dh), h0(dh), G(dh,dh),
     1                 A(dh,dq), B(dh,dv)
      double precision W(dv,dv), GB(dh,dh), 
     1                 m(dh), dm(dh),
     1                 temp51(dh),
     1                 temp31(dv), tempp31(dv),
     1                 chi2ar(1,1)

      integer fvSVDinv
      external fvSVDinv
      integer status 
      integer fvCalcW
      external fvCalcW

      fvRetlif = NORMAL

C -- W
      status = fvCalcW(W,G,B)
      if (status .eq. INFO) write(plun,'(1x,a)')
     1      'editing of singular values done during inv. of W'
      if (iand(status,1).NE.1) then
        fvRetlif = status
        return
      end if
C -- GB
      call fvCalcGB(GB,W,B,G)
C -- Gvp
      call fvsATBA(Gvp, A,GB,dh,dv)
      call fvAMB(Gvp,Gv,Gvp,dv,dv)
C -- Cp = Gvp^(-1)
      call fvCopy(Cp,Gvp,dv,dv)
      status = fvSVDinv(Cp, dv)
      if (status .eq. INFO) write(plun,'(1x,a)')
     1    'editing of singular values done during inv. of C^(-1)'
      if (iand(status,1).NE.1) then
        fvRetlif = status
        return
      end if
C -- m = p0 - h0
      call fvAMB(m,p0,h0,dh,1)
C -- vp = Cp.(Gv.v - A^T.GB.(p0-h0))
      call fvAB(temp31,Gv,v,dv,dv,1)
      call fvAB(temp51,GB,m,dh,dh,1)
      call fvATB(tempp31,A,temp51,dh,dv,1)
      call fvAMB(temp31,temp31,tempp31,dv,1)
      call fvAB(vp,Cp,temp31,dv,dv,1)

C -- chi2
C
C -- chi2 = r^T.G.r + (v-vp)^T.Gvp.(v-vp)
C -- r = p0 - h0 - A.v - B.q
C -- BUT:
C -- shoudn't this be h(vp,qp) at new vertex vp?
C -- i.e. should q be calculated using fvq(q, vp,p0)?
C -- do it...
CCC      call fvq(qp, vp,p0)
C -- Fruehwirt CERN 90-06 says:
C -- distance of track from the new vertex is expressed by the chi-square
C -- of the smoothed residuals
C -- chi2 = (p0-h0-A.v-B.q)^T.G.(p0-h0-A.v-B.q) + (v-vp)^T.Gvp.(v-vp)
      call fvAB(temp51, A,v,dh,dv,1)
      call fvAMB(dm, m,temp51,dh,1)
      call fvAB(temp51, B,q,dh,dv,1)
      call fvAMB(temp51, dm,temp51, dh, 1)
      call fvsATBA(chi2ar, temp51,G,dh,1)
      chi2 = chi2ar(1,1)

      call fvAMB(temp31, v, vp, dv, 1)
      call fvsATBA(chi2ar, temp31,Gvp,dv,1)
      chi2 = chi2 + chi2ar(1,1)

      return
      end
      integer function fvv0q0(v0, Cv0, q0, D0, x0,C x0, p0, Cp0)
C--------------------------------------------------------------------------
C -- calculate and return v0, Cv0, q0, D0 from x0, Cx0, p0, Cp0
C -- 
C -- calculate vertex position v0 and 3-momentum vector q0 
C -- from track parameters p0 (reference surface: r = d0)
C -- at carthesian coordinates x0
C -- 
C -- p0 = {w, tl, psi0, d0, z0}, x0 = {x, y, z}
C --
C -- ????
C -- we do not actually calculate D0, as it is not used
C -- otherwise it should be 1/delta, delta -> inf.
C -- ????
C --
C
      implicit none
      include 'fv.inc'
      double precision v0(dv), Cv0(dv,dv), q0(dq), D0(dq,dq),
     1                 x0(dv), Cx0(dv,dv), p0(dh), Cp0(dh,dh)
      integer status, fvq, fvvr, fvCvr
      external fvq, fvvr, fvCvr

      fvv0q0 = NORMAL
C -- calculate vertex vector in appropriate coord. system
      status = fvvr(v0, x0)
      if (iand(status,1).NE.1) then  
        fvv0q0 = status
        return
      end if
C -- calculate direction vector q0 at min dist. of helix p0 to point v0
      status = fvq(q0, v0,p0)
      if (iand(status,1).NE.1) then  
        fvv0q0 = status
        return
      end if

C -- tranform covariance matrix
      status = fvCvr(Cv0, v0,x0,Cx0)
      if (iand(status,1).NE.1) then  
        fvv0q0 = status
        return
      end if

      return
      end
      integer function fvCvr(Cv, v,x,Cx)
C----------------------------------------------------------------------
C -- transform covariance matrix Cx of x to get Cv
C --
      implicit none
      include 'fv.inc'
      double precision Cv(dv,dv), v(dv),x(dv),Cx(dv,dv)

      double precision r, phi, z, J(dv,dv)

      fvCvr = NORMAL
      r = v(1)
      phi = v(2)
      z = v(3)

C -- Jacobi
C -- J = d(r,phi,z)/d(x,y,z)
      J(1, 1) = (x(1)-RR0)/r
      J(1, 2) = x(2)/r
      J(1, 3) = 0d0

      J(2, 1) = -x(2)/r/r
      J(2, 2) = (x(1)-RR0)/r/r
      J(2, 3) = 0d0

      J(3, 1) = 0d0
      J(3, 2) = 0d0
      J(3, 3) = 1d0

C -- Cv = J.Cx.J^T
      call fvsABAT(Cv,J,Cx,dv,dv)

      return
      end
      integer function fvvr(v, x)
C--------------------------------------------------------------------------
C -- calculate vertex vector v = {r,phi,z} from x = {x, y, z}
C -- i.e. transform to a cylindrical coordinate system
C --
C
      implicit none
      include 'fv.inc'
      double precision v(dv), x(dv)
      double precision xx,yy
      
      fvvr = NORMAL
      xx = x(1)-RR0
      yy = x(2)
      v(1) = sqrt(xx*xx + yy*yy)
      if (v(1) .eq. 0d0) then  
C -- can't handle vertex at (0,0,z) at the moment!
        fvvr = ERROR
        return
      end if
      v(2) = atan2(yy,xx)
      v(3) = x(3)

      return
      end 
      integer function fvx(x, v)
C--------------------------------------------------------------------------
C -- calculate vector x = {x, y, z} from v = {r,phi,z}
C -- i.e. transform back from cylindrical coordinate system
C --
C
      implicit none
      include 'fv.inc'
      double precision x(dv), v(dv)
      double precision r, phi, z, sphi, cphi
      
      fvx = NORMAL

      r = v(1)
      phi = v(2)
      z = v(3)

      sphi = sin(phi)
      cphi = cos(phi)
      x(1) = r*cphi+RR0
      x(2) = r*sphi
      x(3) = z

      return
      end 
      integer function fvxp(x,Cx,p,Cp,v,q,C,D,E,A,B)
C--------------------------------------------------------------------------
C -- calculate and return x, Cx, p, Cp from v, q, C, D, E, A, B
C -- A and B given for conveniance only
C
      implicit none
      include 'fv.inc'
      double precision x(dv), Cx(dv,dv), p(dh), Cp(dh,dh),
     1                 v(dv), q(dv), C(dv,dv), D(dq,dq), E(dv,dq),
     1                 A(dh,dv), B(dh,dq)
      double precision r, phi, z, w, tl, psi,
     1                 oow, xi, cxi, sxi, gamma,
     1                 sphi, cphi, J(dv,dv), temp55(dh,dh)

      fvxp = NORMAL

      r = v(1)
      phi = v(2)
      z = v(3)
      
      w = q(1)
      tl = q(2)
      psi = q(3)

      oow = 1d0/w
      xi = psi - phi
      cxi = cos(xi)
      sxi = sin(xi)
      gamma = atan(r*cxi/(oow-r*sxi))

      sphi = sin(phi)
      cphi = cos(phi)
      x(1) = r*cphi+RR0
      x(2) = r*sphi
      x(3) = z

      p(1) = w
      p(2) = tl
      p(3) = psi - gamma
      p(4) = -r*cxi/sin(gamma)+oow
      p(5) = z - gamma/w*tl

C -- transform covariance matrix
C -- J = d(x,y,z)/d(r,phi,z)
      J(1, 1) = cphi
      J(1, 2) = -r*sphi
      J(1, 3) = 0d0

      J(2, 1) = sphi
      J(2, 2) = r*cphi
      J(2, 3) = 0d0

      J(3, 1) = 0d0
      J(3, 2) = 0d0
      J(3, 3) = 1d0

C -- Cx = J.C.J^T
      call fvsABAT(Cx,J,C,dv,dv)

C -- Cp = A.C.A^T + B.E.A^T + A.E^T.B^T + B.D.B^T
      call fvsABAT(Cp,A,C,dh,dv)
      call fvABCT(temp55,B,E,A,dh,dv)
      call fvAPB(Cp,Cp,temp55,dh,dh)
      call fvABTCT(temp55,A,E,B,dh,dv)
      call fvAPB(Cp,Cp,temp55,dh,dh)
      call fvsABAT(temp55,B,D,dh,dv)
      call fvAPB(Cp,Cp,temp55,dh,dh)

      return
      end
      integer function fvq(q, v,h)
C--------------------------------------------------------------------------
C -- calculate 3-momentum vector q = {w, tl, psi}
C -- from helix parameters h (reference surface: r = d0)
C -- at minimum distance to space point v
C -- 
C -- h = {w, tl, psi0, d0, z0}, v = {x, y, z}
C --
C
      implicit none
      include 'fv.inc'
      double precision q(dq), v(dv), h(dh)
      double precision w, tl, psi0, gamma, psi
CCC      double precision dd0, r0, u0, xi, gamma, psi, r, phi 

      fvq = NORMAL

CCC      r = sqrt(v(1)*v(1) + v(2)*v(2))
CCC      phi =  atan2(v(2),v(1))

      w =     h(1)
      tl =    h(2)
      psi0 =  h(3)
CCC      dd0 =   h(4)
C
C -- calculate psi at a position on track helix p, where
C -- the distance between point v and the helix is minimal 
C -- in the x-y plane. 
C -- How can we include the distance in z? Is this necessary?
C
C -- for the moment, just return track parameters in global system,
C -- i.e. at minimum distance of approach to origin
CCC      if (dd0 .EQ. 0d0) then
        gamma = 0d0
CCC      else 
CCC        r0 = abs(dd0)
CCC        u0 = 2d0*r0*atan(-dd0*cos(psi0)/(r0+dd0*sin(psi0)))
CCC        xi = -u0/dd0 - phi
CCC        gamma = atan(r*sin(xi)/(1d0/w-dd0 - r*cos(xi)))
CCC      end if

      psi = psi0 + gamma

      q(1) = w
      q(2) = tl
      q(3) = psi

      return
      end
      integer function fvrh2q(q,Cq, r,h,Ch)
C--------------------------------------------------------------------------
C -- calculate 3-momentum vector q = {w, tl, psi} and covariance matrix Cq
C -- from track parameters h = {w, tl, psi0, d0, z0} (reference surface is
C -- Rref = d0) and covariance matrix Ch at a given radius r > d0
C --
C
      implicit none
      include 'fv.inc'
      double precision q(dq),Cq(dq,dq), r,h(dh),Ch(dh,dh)
      double precision w, tl, psi0, d0,
     1                 gamma, psi, sdtw,srmd,dgamdu,J(dq,dh)
      integer i0,i1

      fvrh2q = NORMAL

      w =     h(1)
      tl =    h(2)
      psi0 =  h(3)
      d0 =   h(4)

C -- track crosses r only, if
C -- d0<r and d0*w < 1
      if (abs(d0) .LT. r) then
        srmd = sqrt(r*r-d0*d0)
        sdtw = sqrt(1-d0*w)

        gamma = 2d0*asin(w/2d0*srmd/sdtw)
        dgamdu = 0.5d0/srmd/sdtw**3/sqrt(1d0-(srmd*w/2d0/sdtw)**2)

        psi = psi0 + gamma

        q(1) = w
        q(2) = tl
        q(3) = psi

        J(1,1) = 1d0
        J(1,2) = 0d0
        J(1,3) = 0d0
        J(1,4) = 0d0
        J(1,5) = 0d0

        J(2,1) = 0d0
        J(2,2) = 1d0
        J(2,3) = 0d0
        J(2,4) = 0d0
        J(2,5) = 0d0

        J(3,1) = (2d0 - d0*w)*srmd*srmd*dgamdu
        J(3,2) = 0d0
        J(3,3) = 1d0
        J(3,4) = w*(-2d0*d0 + d0*d0*w + r*r*w)*dgamdu
        J(3,5) = 1d0

        call fvsABAT(Cq, J,Ch,dq,dh)
      else
        psi = psi0
        do i0=1,dq
        do i1=1,dq
          Cq(i0,i1) = Ch(i0,i1)
        end do
        end do
      end if

      return
      end
      integer function fvrh2v(v,Cv, r,h,Ch)
C--------------------------------------------------------------------------
C -- calculate vertex vector v = {vx, vy, vz} and covariance matrix Cv
C -- from track parameters h = {w, tl, psi0, d0, z0} (reference surface is
C -- Rref = d0) and covariance matrix Ch at a given radius r > d0
C
C---??????????????????????not yet done
C
C --
C
      implicit none
      include 'fv.inc'
      double precision v(dv),Cv(dv,dv), r,h(dh),Ch(dh,dh)
      double precision w, tl, psi0, d0, z0, xx,yy,zz,
     1                 gamma, psi, sdtw,srmd,dgamdu,J(dq,dh)
      integer i0,i1

      fvrh2v = NORMAL

      w =     h(1)
      tl =    h(2)
      psi0 =  h(3)
      d0 =    h(4)
      z0 =    h(5)

C -- track crosses r only, if
C -- d0<r and d0*w < 1
      if (abs(d0) .LT. r) then
        srmd = sqrt(r*r-d0*d0)
        sdtw = sqrt(1-d0*w)

        gamma = 2d0*asin(w/2d0*srmd/sdtw)
        dgamdu = 0.5d0/srmd/sdtw**3/sqrt(1d0-(srmd*w/2d0/sdtw)**2)
        psi = psi0 + gamma

        v(1) = xx
        v(2) = yy
        v(3) = zz

        J(1,1) = 1d0
        J(1,2) = 0d0
        J(1,3) = 0d0
        J(1,4) = 0d0
        J(1,5) = 0d0

        J(2,1) = 0d0
        J(2,2) = 1d0
        J(2,3) = 0d0
        J(2,4) = 0d0
        J(2,5) = 0d0

        J(3,1) = (2d0 - d0*w)*srmd/
     1           (2d0*sdtw**3*sqrt(1d0-(srmd*w/2d0/sdtw)**2))
        J(3,2) = 0d0
        J(3,3) = 1d0
        J(3,4) = w*(-2d0*d0 + d0*d0*w + r*r*w)/
     1           (2d0*srmd*sdtw**3*
     1            sqrt(1d0-(srmd*w/2d0/sdtw)**2))
        J(3,5) = 1d0

        call fvsABAT(Cv, J,Ch,dv,dh)
      else
        psi = psi0
        do i0=1,dv
        do i1=1,dv
          Cv(i0,i1) = Ch(i0,i1)
        end do
        end do
      end if

      return
      end
      integer function fvhCh(h,Ch,v,q,C,D,E)
C--------------------------------------------------------------------------
C -- calculate and return helix parameters h, Ch from v, q, C, D, E
C
      implicit none
      include 'fv.inc'
      double precision h(dh), Ch(dh,dh),
     1                 v(dv), q(dq), C(dv,dv), D(dq,dq), E(dv,dq)

      integer fvCh, fvh, status
      external fvCh, fvh

      fvhCh = NORMAL
      status = fvh(h, v,q)
      if (iand(status,1).NE.1) fvhCh = status
      
      status = fvCh(Ch, v,q,C,D,E)
      if (iand(status,1).NE.1) fvhCh = status

      return
      end
      integer function fvh(h, v,q)
C--------------------------------------------------------------------------
C -- calculate and return helix parameters h={w,tl,psi0,d0,z0} 
C -- from v={vx,vy,vz} and q={w,tl,psi}
C
      implicit none
      include 'fv.inc'
      double precision h(dh),
     1                 v(dv), q(dq)
      double precision xx,yy,r, phi, z, w, tl, psi,
     1                 oow, xi, cxi, sxi, gamma

      fvh = NORMAL

      xx = v(1)
      yy = v(2)

      r = sqrt(xx*xx+yy*yy)
      phi =  atan2(yy,xx)
      z = v(3)
      
      w = q(1)
      tl = q(2)
      psi = q(3)

      if (w .NE. 0d0) then
        oow = 1d0/w
        xi = psi - phi
        cxi = cos(xi)
        sxi = sin(xi)
        gamma = atan(r*cxi/(oow-r*sxi))
        h(1) = w
        h(2) = tl
        h(3) = psi - gamma
        h(4) = oow - (oow-r*sxi)/cos(gamma)
        h(5) = z - gamma/w*tl
      else
C -- for neutral particles
C -- ?????? check this out
        xi = psi - phi
        sxi = sin(xi)
        h(1) = w
        h(2) = tl
        h(3) = psi
        h(4) = r*sxi
        h(5) = z
      end if

      write(plun,'(1x,a,5(g10.3))') "fvh",h
      write(plun,'(1x,a,3(g10.3))') "v",v
      write(plun,'(1x,a,3(g10.3))') "q",q
      write(plun,'(1x,a,g10.3)') "gamma",gamma
      return
      end
      integer function fvCh(Ch, v,q,C,D,E)
C--------------------------------------------------------------------------
C -- calculate and return cavariance matrix of helix parameters 
C -- Ch = cov{Ch} from covariances of vertex v={vx,vy,vz} and momentum 
C -- q={w,tl,psi},
C -- C = cov{v}, D = cov{q}, E=cov{v,q}
C
      implicit none
      include 'fv.inc'
      double precision Ch(dh,dh),
     1                 v(dv),q(dq),C(dv,dv), D(dq,dq), E(dv,dq)
      double precision temp55(dh,dh),
     1                 A(dh,dv), B(dh,dq), h0(dh)
      integer fvABh0
      external fvABh0

C -- Ch = A.C.A^T + B.E.A^T + A.E^T.B^T + B.D.B^T
      fvCh = fvABh0(A,B,h0, v,q)

      call fvsABAT(Ch,A,C,dh,dv)
      call fvABCT(temp55,B,E,A,dh,dv)
      call fvAPB(Ch,Ch,temp55,dh,dh)
      call fvABTCT(temp55,A,E,B,dh,dv)
      call fvAPB(Ch,Ch,temp55,dh,dh)
      call fvsABAT(temp55,B,D,dh,dv)
      call fvAPB(Ch,Ch,temp55,dh,dh)

      return
      end
      integer function fvABh0(A, B, h0, v0, q0)
C--------------------------------------------------------------------------
C -- calculate and return A, B, h0 from v0, q0:
C -- 
C -- calculate coefficients of taylor expansion for estimated
C -- track parameters p around point (v0,q0):
C --
C --   p    =    h(v,q) + eps
C --     \approx h(v0,q0) + A.(v-v0) + B.(q-q0) + eps
C --        =    h0 + A.v + B.q + eps
C -- where
C --   A  =  dp/dv(v0,q0),
C --   B  =  dp/dq(v0,q0), and 
C --   h0 =  h(v0,q0) - A.v0 - B.q0 are the
C -- derivatives of track parameters w/r to vertex position
C -- and track 3-momentum at space point v0 and estimated momentum q0
C --
C
      implicit none
      include 'fv.inc'
      double precision A(dh,dv), B(dh,dq), h0(dh), 
     1                 v0(dv), q0(dq)
      double precision xx, yy, r, phi, z, w, tl, psi,
     1                 xi, cxi, sxi, oow, rw, gamma, sg, cg,
     1                 psi0, d0, z0,
     1                 drdx,drdy,rdxidx,rdxidy,
     1                 dgdx, dgdy, dgdw, dgdpsi, dgdvar0
      double precision twopi

      parameter (twopi = 6.28318530717959)

      fvABh0 = NORMAL

      xx = v0(1)
      yy = v0(2)
      r = sqrt(xx*xx+yy*yy)
      phi = atan2(yy,xx)
      z = v0(3)

      w = q0(1)
      tl = q0(2)
      psi = q0(3)

C -- calculate some derived quantities
      xi = mod(psi-phi,twopi)
      cxi = cos(xi)
      sxi = sin(xi)
      oow = 1/w
      rw = r*w

      gamma = atan(r*cxi/(oow-r*sxi))
      sg = sin(gamma)
      cg = cos(gamma)

C -- calc transformed quantities
      psi0 = psi - gamma
      d0 = oow - (oow-r*sxi)/cg
      z0 = z - gamma/w*tl

C -- calc Jacobian 

      if (r .NE. 0d0) then
        drdx = xx/r
        drdy = yy/r
        rdxidx = yy/r
        rdxidy = -xx/r
      else
        drdx = 0d0
        drdy = 0d0
        rdxidx = 0d0
        rdxidy = 0d0
      end if

      dgdvar0 = 1d0/(1d0+rw*rw-2d0*rw*sxi)
      dgdx = dgdvar0*(w*cxi*drdx+w*(rw-sxi)*rdxidx)
      dgdy = dgdvar0*(w*cxi*drdy+w*(rw-sxi)*rdxidy)
      dgdw = dgdvar0*r*cxi
      dgdpsi = dgdvar0*rw*(rw-sxi)

C -- d w / d r, d phi, d z
      A(1,1) = 0d0
      A(1,2) = 0d0
      A(1,3) = 0d0
C -- d tl / d x, d y, d z
      A(2,1) = 0d0
      A(2,2) = 0d0
      A(2,3) = 0d0
C -- d psi0 / d x, d y, d z
      A(3,1) = -dgdx
      A(3,2) = -dgdy
      A(3,3) = 0d0
C -- d d0 / d x, d y, d z
      A(4,1) = cxi/cg*rdxidx + sxi/cg*drdx-
     1          (oow-r*sxi)*sg/cg/cg*dgdx
      A(4,2) = cxi/cg*rdxidy + sxi/cg*drdy-
     1          (oow-r*sxi)*sg/cg/cg*dgdy
      A(4,3) = 0d0
C -- d z0 / d x, d y, d z
      A(5,1) = -tl/w*dgdx
      A(5,2) = -tl/w*dgdy
      A(5,3) = 1d0

C -- d w / d w, d tl, d psi
      B(1,1) = 1d0
      B(1,2) = 0d0
      B(1,3) = 0d0
C -- d tl / d w, d tl, d psi
      B(2,1) = 0d0
      B(2,2) = 1d0
      B(2,3) = 0d0
C -- d psi0 / d w, d tl, d psi
      B(3,1) = -dgdw
      B(3,2) = 0d0
      B(3,3) = 1d0-dgdpsi
C -- d d0 / d w, d tl, d psi
      B(4,1) = -oow*oow*(1d0-1d0/cg)-
     1          (oow-r*sxi)/cg/cg*sg*dgdw
      B(4,2) = 0d0
      B(4,3) = r*cxi/cg - 
     1          (oow-r*sxi)/cg/cg*sg*dgdpsi
C -- d z0 / d w, d tl, d psi
      B(5,1) = -tl/w*(dgdw - gamma/w)
      B(5,2) = -gamma/w
      B(5,3) = -tl/w*dgdpsi

      h0(1) = 0d0
      h0(2) = 0d0
      h0(3) = psi0 - A(3,1)*v0(1)-A(3,2)*v0(2)-
     1        B(3,1)*q0(1)-B(3,3)*q0(3)
      h0(4) = d0 - A(4,1)*v0(1)-A(4,2)*v0(2)-
     1        B(4,1)*q0(1)-B(4,3)*q0(3)
      h0(5) = z0 - A(5,1)*v0(1)-A(5,2)*v0(2)-A(5,3)*v0(3 )-
     1        B(5,1)*q0(1)-B(5,2)*q0(2)-B(5,3)*q0(3)

      return
      end
      integer function fvCalcG(G,Cov,n)
C-----------------------------------------------------------------------
C -- calculate Gk = cov(pk)^(-1) for n x n matrix Cov
C
      implicit none
      include 'fv.inc'
      integer n
      double precision G(n,n), Cov(n,n)
      integer fvSVDinv
      external fvSVDinv

C -- invert covariance matrix
      call fvCopy(G,Cov,n,n)
      fvCalcG = fvSVDinv(G,n)

      return
      end
      integer function fvCalcW(W,G,B)
C-----------------------------------------------------------------------
C -- calculate W = (B^T.G.B)^(-1)
C
      implicit none
      include 'fv.inc'
      double precision W(dv,dv), B(dh,dv), G(dh,dh)
      integer fvSVDinv
      external fvSVDinv

      call fvsATBA(W, B,G,dh,dv)
      fvCalcW = fvSVDinv(W, dv)

      return
      end
      subroutine fvCalcGB(GB,W,B,G)
C-----------------------------------------------------------------------
C -- calculate GB = G - G.B.W.B^T.G
C
      implicit none
      include 'fv.inc'
      double precision GB(dh,dh),W(dv,dv), B(dh,dv), G(dh,dh)

      double precision temp55(dh,dh)

      call fvsABAT(temp55,B,W,dh,dv)
      call fvsABAT(GB,G,temp55,dh,dh)
      call fvAMB(GB,G,GB,dh,dh)

      return
      end
      subroutine fvError(text)
C-----------------------------------------------------------------------
C -- output error string
C
      implicit none
      character *(*) text
      write(*,'(1x,a)') text
      return
      end
      subroutine fvSPrt(p,l)
C-----------------------------------------------------------------------
C -- set print flag and print unit
C
      implicit none
      include 'fv.inc'
      logical p
      integer l
      print = p
      plun = l
      return
      end
      subroutine fvHelix2P4(t,Ct, w2pt,h,Ch)
C----------------------------------------------------------------------- 
C -- calculate 3-momentum from helix parameters parameters
C -- fill into 4-momentum array
C -- use w2pt for calculating pt from curvature w
C
      implicit none
      include 'fv.inc'

      double precision t(4), Ct(4,4)
      double precision w2pt,h(dh), Ch(dh,dh)


      double precision w, tl, psi0
      double precision pt, sph, cph, px, py, pz
      double precision ps, dpdk, xy,
     1  sxx, sxy, syy, sxz, syz, szz

      w = h(1)
      tl = h(2)
      psi0 = h(3)

      sph = sin(psi0)
      cph = cos(psi0)
      pt = w2pt/abs(w)

      px = pt*cph
      py = pt*sph
      pz = pt*tl

C -- calculate track momentum 3-vector and its error matrix
      ps = w2pt / w
      dpdk = ps*ps/w2pt
      xy = 2.*ps*dpdk*cph*sph*Ch(1,3)

      sxx = (dpdk*cph)**2 * Ch(1,1) + (ps*sph)**2 * Ch(3,3) + xy
      sxy = cph*sph*(dpdk*dpdk*Ch(1,1) - ps*ps*Ch(3,3)) +
     1      ps*dpdk*(sph*sph-cph*cph)*Ch(1,3)
      syy = (dpdk*sph)**2 * Ch(1,1) + (ps*cph)**2 * Ch(3,3) - xy
      sxz = dpdk*dpdk*cph*tl*Ch(1,1) -
     1      ps*dpdk*(cph*Ch(1,2)-sph*tl*Ch(1,3)) -
     1      ps*ps*sph*Ch(2,3)
      syz = dpdk*dpdk*sph*tl*Ch(1,1) -
     1      ps*dpdk*(sph*Ch(1,2) + cph*tl*Ch(1,3)) +
     1      ps*ps*cph*Ch(2,3)
      szz = (dpdk*tl)**2 * Ch(1,1) + ps*ps*Ch(2,2) - 
     1      2.*ps*dpdk*tl*Ch(1,2)

      t(1) = px
      t(2) = py
      t(3) = pz

      Ct(1,1) = sxx
      Ct(1,2) = sxy
      Ct(1,3) = sxz

      Ct(2,1) = sxy
      Ct(2,2) = syy
      Ct(2,3) = syz

      Ct(3,1) = sxz
      Ct(3,2) = syz
      Ct(3,3) = szz

      return
      end
      subroutine fvHelix2P3(t,Ct, w2pt,h,Ch)
C----------------------------------------------------------------------- 
C -- calculate 3-momentum from helix parameters h, Ch
C -- fill into 3-momentum array
C -- use w2pt for calculating pt from curvature w
C
      implicit none
      include 'fv.inc'

      double precision t(dv), Ct(dv,dq)
      double precision w2pt,h(dh),Ch(dh,dh)


      double precision w, tl, psi0
      double precision pt, sph, cph, px, py, pz
      double precision ps, dpdk, xy,
     1  sxx, sxy, syy, sxz, syz, szz

      w = h(1)
      tl = h(2)
      psi0 = h(3)

      sph = sin(psi0)
      cph = cos(psi0)
      pt = w2pt/abs(w)

      px = pt*cph
      py = pt*sph
      pz = pt*tl

C -- calculate track momentum 3-vector and its error matrix
      ps = w2pt / w
      dpdk = ps*ps/w2pt
      xy = 2.*ps*dpdk*cph*sph*Ch(1,3)

      sxx = (dpdk*cph)**2 * Ch(1,1) + (ps*sph)**2 * Ch(3,3) + xy
      sxy = cph*sph*(dpdk*dpdk*Ch(1,1) - ps*ps*Ch(3,3)) +
     1      ps*dpdk*(sph*sph-cph*cph)*Ch(1,3)
      syy = (dpdk*sph)**2 * Ch(1,1) + (ps*cph)**2 * Ch(3,3) - xy
      sxz = dpdk*dpdk*cph*tl*Ch(1,1) -
     1      ps*dpdk*(cph*Ch(1,2)-sph*tl*Ch(1,3)) -
     1      ps*ps*sph*Ch(2,3)
      syz = dpdk*dpdk*sph*tl*Ch(1,1) -
     1      ps*dpdk*(sph*Ch(1,2) + cph*tl*Ch(1,3)) +
     1      ps*ps*cph*Ch(2,3)
      szz = (dpdk*tl)**2 * Ch(1,1) + ps*ps*Ch(2,2) - 
     1      2.*ps*dpdk*tl*Ch(1,2)

      t(1) = px
      t(2) = py
      t(3) = pz

      Ct(1,1) = sxx
      Ct(1,2) = sxy
      Ct(1,3) = sxz

      Ct(2,1) = sxy
      Ct(2,2) = syy
      Ct(2,3) = syz

      Ct(3,1) = sxz
      Ct(3,2) = syz
      Ct(3,3) = szz

      return
      end
      subroutine fvQ2P4(p,Cp, q,Cq,w2pt)
C--------------------------------------------------------------------------
C -- calculate and return 3-momentum in four-vector p 
C -- and its cov. matrix Cp in carth. coord. from 
C -- q = {w, tl, psi}
C -- note: p, Cp is 4x4!
C --
C
      implicit none
      include 'fv.inc'
      double precision p(4),Cp(4,4),q(dq),Cq(dq,dq),w2pt


      double precision w, tl, p0
      double precision Cp11, Cp12, Cp13, Cp22, Cp23, Cp33

      double precision pt, sph, cph, px, py, pz
      double precision ps, dpdk, xy,
     1  sxx, sxy, syy, sxz, syz, szz

      w = q(1)
      tl = q(2)
      p0 = q(3)

      sph = sin(p0)
      cph = cos(p0)
      pt = w2pt/abs(w)

      Cp11 = Cq(1,1)
      Cp12 = Cq(1,2)
      Cp13 = Cq(1,3)
      Cp22 = Cq(2,2)
      Cp23 = Cq(2,3)
      Cp33 = Cq(3,3)

      px = pt*cph
      py = pt*sph
      pz = pt*tl

C -- calculate track momentum 3-vector and its error matrix
      ps = w2pt / w
      dpdk = ps*ps/w2pt
      xy = 2.*ps*dpdk*cph*sph*Cp13

      sxx = (dpdk*cph)**2 * Cp11 + (ps*sph)**2 * Cp33 + xy
      sxy = cph*sph*(dpdk*dpdk*Cp11 - ps*ps*Cp33) +
     1      ps*dpdk*(sph*sph-cph*cph)*Cp13
      syy = (dpdk*sph)**2 * Cp11 + (ps*cph)**2 * Cp33 - xy
      sxz = dpdk*dpdk*cph*tl*Cp11 -
     1      ps*dpdk*(cph*Cp12-sph*tl*Cp13) -
     1      ps*ps*sph*Cp23
      syz = dpdk*dpdk*sph*tl*Cp11 -
     1      ps*dpdk*(sph*Cp12 + cph*tl*Cp13) +
     1      ps*ps*cph*Cp23
      szz = (dpdk*tl)**2 * Cp11 + ps*ps*Cp22 - 
     1      2.*ps*dpdk*tl*Cp12

      p(1) = px
      p(2) = py
      p(3) = pz

      Cp(1,1) = sxx
      Cp(1,2) = sxy
      Cp(1,3) = sxz

      Cp(2,1) = sxy
      Cp(2,2) = syy
      Cp(2,3) = syz

      Cp(3,1) = sxz
      Cp(3,2) = syz
      Cp(3,3) = szz

      return
      end
      subroutine fvPE(p,Cp, it,w2pt,m,ql,Cql)
C-----------------------------------------------------------------------        
C -- calculate 4-momentum for track it
      implicit none
      include 'fv.inc'
      double precision p(4),Cp(4,4)
      integer it
      double precision w2pt,m,ql(dq,*),Cql(dq,dq,*)

      call fvQ2P4(p,Cp, ql(1,it),Cql(1,1,it),w2pt)
      call fvEnergy(p, Cp, m)

      return
      end
      subroutine fvEnergy(t,Ct,m)
C-----------------------------------------------------------------------        
C -- calculate energy from 3-momentum
C -- fill entries in t and Ct
C
      implicit none
      include 'fv.inc'

      double precision t(4), Ct(4,4), m

      double precision px,py,pz,e

      px = t(1)
      py = t(2)
      pz = t(3)
      e = sqrt(m*m+ px*px+py*py+pz*pz)
      t(4) = e

      Ct(1,4) = (px*Ct(1,1) + py*Ct(1,2) + pz*Ct(1,3))/e
      Ct(2,4) = (px*Ct(1,2) + py*Ct(2,2) + pz*Ct(2,3))/e
      Ct(3,4) = (px*Ct(1,3) + py*Ct(2,3) + pz*Ct(3,3))/e
      Ct(4,4) = (px*px*Ct(1,1) + py*py*Ct(2,2) + pz*pz*Ct(3,3) + 
     1      2.*(px*(py*Ct(1,2) + pz*Ct(1,3)) + py*pz*Ct(2,3)))/e/e

      Ct(4,1) = Ct(1,4)
      Ct(4,2) = Ct(2,4)
      Ct(4,3) = Ct(3,4)

      return
      end
      subroutine fvQ2P3(p,Cp, q,Cq,w2pt)
C--------------------------------------------------------------------------
C -- calculate and return 3-momentum vector p
C -- and its cov. matrix Cp in carth. coord. from 
C -- q = {w, tl, psi}
C -- note: p, Cp is 3x3!
C --
C
      implicit none
      include 'fv.inc'
      double precision p(dq),Cp(dq,dq),q(dq),Cq(dq,dq),w2pt


      double precision w, tl, p0
      double precision Cp11, Cp12, Cp13, Cp22, Cp23, Cp33

      double precision pt, sph, cph, px, py, pz
      double precision ps, dpdk, xy,
     1  sxx, sxy, syy, sxz, syz, szz

      w = q(1)
      tl = q(2)
      p0 = q(3)

      sph = sin(p0)
      cph = cos(p0)
      pt = w2pt/abs(w)

      Cp11 = Cq(1,1)
      Cp12 = Cq(1,2)
      Cp13 = Cq(1,3)
      Cp22 = Cq(2,2)
      Cp23 = Cq(2,3)
      Cp33 = Cq(3,3)

      px = pt*cph
      py = pt*sph
      pz = pt*tl

C -- calculate track momentum 3-vector and its error matrix
      ps = w2pt / w
      dpdk = ps*ps/w2pt
      xy = 2.*ps*dpdk*cph*sph*Cp13

      sxx = (dpdk*cph)**2 * Cp11 + (ps*sph)**2 * Cp33 + xy
      sxy = cph*sph*(dpdk*dpdk*Cp11 - ps*ps*Cp33) +
     1      ps*dpdk*(sph*sph-cph*cph)*Cp13
      syy = (dpdk*sph)**2 * Cp11 + (ps*cph)**2 * Cp33 - xy
      sxz = dpdk*dpdk*cph*tl*Cp11 -
     1      ps*dpdk*(cph*Cp12-sph*tl*Cp13) -
     1      ps*ps*sph*Cp23
      syz = dpdk*dpdk*sph*tl*Cp11 -
     1      ps*dpdk*(sph*Cp12 + cph*tl*Cp13) +
     1      ps*ps*cph*Cp23
      szz = (dpdk*tl)**2 * Cp11 + ps*ps*Cp22 - 
     1      2.*ps*dpdk*tl*Cp12

      p(1) = px
      p(2) = py
      p(3) = pz

      Cp(1,1) = sxx
      Cp(1,2) = sxy
      Cp(1,3) = sxz

      Cp(2,1) = sxy
      Cp(2,2) = syy
      Cp(2,3) = syz

      Cp(3,1) = sxz
      Cp(3,2) = syz
      Cp(3,3) = szz

      return
      end
      subroutine fvp3(p, q,w2pt)
C--------------------------------------------------------------------------
C -- calculate and return 3-momentum vector p
C -- Cp in carth. coord. from 
C -- q = {w, tl, psi}
C --
C
      implicit none
      include 'fv.inc'
      double precision p(dq),q(dq),w2pt


      double precision w, tl, p0
      double precision pt, sph, cph

      w = q(1)
      tl = q(2)
      p0 = q(3)

      sph = sin(p0)
      cph = cos(p0)
      pt = w2pt/abs(w)

      p(1) = pt*cph
      p(2) = pt*sph
      p(3) = pt*tl

      return
      end
      subroutine fvq2pvec(p,q,w2pt)
C--------------------------------------------------------------------------
C -- calculate and return 3-momentum vector p
C -- q = {w, tl, psi}, p = {px,py,pz}
C -- note: p, Cp is 3x3!
C --
C
      implicit none
      include 'fv.inc'
      double precision p(dq),q(dq),w2pt


      double precision pt

      pt = w2pt/abs(q(1))

      p(1) = pt*cos(q(3))
      p(2) = pt*sin(q(3))
      p(3) = pt*q(2)

      return
      end
      subroutine fvp2q(q, p,c,w2pt)
C--------------------------------------------------------------------------
C -- calculate q from 3-momentum vector p
C -- p = {px,py,pz}, q = {w, tl, psi}
C -- c is charge
C --
C
      implicit none
      include 'fv.inc'
      double precision p(dq),q(dq),c,w2pt


      double precision pt

      pt = sqrt(p(1)**2 + p(2)**2)
C -- ??????this is a kludge
      if (c .NE. 0d0) then
        q(1) = c*w2pt/pt
      else
        q(1) = w2pt/pt
      end if
      q(2) = p(3)/pt
      q(3) = atan2(p(2),p(1))

      return
      end
      integer function fvInv6s(C,D,E)
C--------------------------------------------------------------------------
C -- invert symmetric 6 x 6 matrix
C --       / C   E \
C --   M = |  T     |
C --       \ E   D /
C
      implicit none
      include 'fv.inc'
      double precision C(dv,dv), D(dq,dq), E(dv,dq)
      double precision M(dv+dq,dv+dq)
      integer i0,i1
      integer fvSVDinv
      external fvSVDinv

C -- fill up 6x6 matrix
      do i0 = 1, 3
      do i1 = 1, 3
         M(i0,i1) = C(i0,i1)
         M(i0+3,i1+3) = D(i0,i1)
         M(i0,i1+3) = E(i0,i1)
         M(i1+3,i0) = E(i0,i1)
      end do
      end do

      fvInv6s = fvSVDinv(M,dv+dq)

C -- get from 6x6 matric
      do i0 = 1, 3
      do i1 = 1, 3
         C(i0,i1) = M(i0,i1)
         D(i0,i1) = M(i0+3,i1+3)
         E(i0,i1) = M(i0,i1+3)
      end do
      end do

      return
      end
CC      integer function fvInv(A, n)
CCCCC--------------------------------------------------------------------------
CCC -- invert symmetric n x n matrix A of order n, and return A^(-1)
CCC -- A will be overwritten!
CCC --
CCC
CC      implicit none
CC      include 'fv.inc'
CC      integer n
CC      double precision A(n,n)
CC      integer ier
CC      integer temp(20)
CC
CC      fvInv = NORMAL
CCC#ifndef CRAY
CC      call DINV(n, A, n, temp, ier)
CCC      call DSINV(n, A, n, ier)
CCC#else
CCC      call RINV(n, A, n, temp, ier)
CCC#endif
CC      if (ier .ne. 0) then
CC         fvInv = ERROR
CC         return
CC      end if
CC      return
CC      end
      integer function fvSVDinv(a, n)
C--------------------------------------------------------------------------
C -- invert symmetric n x n matrix A of order n, and return A^(-1)
C -- A will be overwritten!
C -- use a double precision version of the routine described in
C -- Press et al., Numerical Recipes, Cambridge Univ. Press
C
      implicit none
      integer n
      double precision a(n,n)
      include 'fv.inc'

      double precision TOL
      parameter (TOL = 1d-14)
      integer NMAX
      parameter (NMAX=20)
      double precision U(NMAX,NMAX), W(NMAX), V(NMAX,NMAX)
      integer i, j, jj
      double precision s,tmp(NMAX, NMAX)
      double precision Jacob(NMAX)

      double precision wMax, thresh

      integer status, fvSVD
      external fvSVD

      logical WITHU, WITHV
      data WITHU/.true./,WITHV/.true./

      fvSVDinv = NORMAL
      if (n .gt. NMAX) then
        fvSVDinv = ERROR
        return
      end if

      do i = 1, n
        Jacob(i) = sqrt(abs(a(i,i)))
      end do

      do i = 1, n
      do j = 1, n
        U(i,j) = a(i,j)/Jacob(i)/Jacob(j)
      end do
      end do

      status = fvSVD(U,n,n,NMAX,NMAX,W,V)
      if (status .ne. NORMAL) then
         fvSVDinv = status
         return
      end if

C -- edit singular values
      wMax = 0.
      do i = 1, n
        if (W(i).gt.wMax) wMax = W(i)
      end do
      thresh = TOL*wMax
      do i = 1, n
	if (W(i) .lt. thresh) then
          write(*,*) 'fvSVDinv: editing singular values'
          write(*,*) '----> ', W(i)
          W(i) = 0d0
          fvSVDinv = INFO
        end if
      end do

C -- compose inverse matrix
C -- A <- A^(-1) = V.W^(-1).U^T

      do i= 1, n
        do j = 1, n
          if (W(i) .NE. 0d0)then
	    tmp(i,j) = U(j,i)/W(i)
          else 
            tmp(i,j) = 0d0
          end if
        end do
      end do

      do i = 1, n
	do j = 1, n
	    s = 0d00
	    do jj = 1, n
		s = s + V(i,jj)*tmp(jj,j)
            end do
	    a(i,j) = s/Jacob(i)/Jacob(j)
        end do
      end do

      return
      end
      integer function fvLUinv(y, a,n)
C--------------------------------------------------------------------------
C -- invert n x n matrix A of order n, and return A^(-1) in Y
C -- using LU decomposition/backsubstitution
C -- A will be overwritten!
C -- use a double precision version of the routine described in
C -- Press et al., Numerical Recipes, Cambridge Univ. Press
C
      implicit none
      integer n
      double precision y(n,n),a(n,n)
      include 'fv.inc'

      integer NMAX
      parameter (NMAX=6)
      double precision indx(NMAX), d,tmp(NMAX,NMAX)
      double precision Jacob(NMAX)
      integer status,i,j,fvLUdcmp
      external fvLUdcmp

      fvLUinv = NORMAL

      if (n.gt.NMAX) then
         fvLUinv = ERROR
         return
      end if

      do i = 1, n
        Jacob(i) = sqrt(abs(a(i,i)))
      end do

      do i=1,n
         do j=1,n
            y(i,j)=0d0
            tmp(i,j) = a(i,j)/Jacob(i)/Jacob(j)
         end do
         y(i,i)=1d0
      end do

C -- decompose the matrix just once
      status = fvLUdcmp(tmp,n,NMAX,indx,d)
      if (iand(status,1).NE.1) then
         fvLUinv = status
         return
      end if
C -- find invers by column
      do j=1,n
        call fvLUbksb(tmp,n,NMAX,indx,y(1,j))
      end do

      do i = 1, n
	do j = 1, n
	    y(i,j) = y(i,j)/Jacob(i)/Jacob(j)
        end do
      end do

      return
      end
      subroutine fvCopy(A, B, n, m)
C--------------------------------------------------------------------------
C -- copy n x m matrix B into A
C
      implicit none
      integer n,m
      double precision A(n,m), B(n,m)
      include 'fv.inc'
      integer i0,i1

C#ifndef CRAY
C      call UCOPY(B, A, n*m*floatSize)
C#else
      do i0 = 1, n
      do i1 = 1, m
        A(i0,i1) = B(i0,i1)
      end do
      end do
C#endif
      return
      end
      subroutine fvZeroA(A, n, m)
C--------------------------------------------------------------------------
C -- zero contents  into A for n x m matrices A
C
      implicit none
      integer n, m
      double precision A(n, m)
      include 'fv.inc'
      integer i, j

      do i = 1, n
      do j = 1, m
        A(i,j) = 0d0
      end do
      end do

      return
      end
      subroutine fvAMB(X, A, B, n, m)
C--------------------------------------------------------------------------
C -- calculate X = A - B for n x m matrices X, A, B
C -- it is allowed for result matrix X to be A or B
C
      implicit none
      integer n, m
      double precision X(n,m), A(n, m), B(n, m)
      include 'fv.inc'
      integer i, j

      do i = 1, n
      do j = 1, m
        X(i,j) = A(i,j) - B(i,j)
      end do
      end do

      return
      end
      subroutine fvNegA(X, A, n, m)
C--------------------------------------------------------------------------
C -- calculate X = -A for n x m matrices X, A
C -- it is allowed for result matrix X to be A
C
      implicit none
      integer n, m
      double precision X(n,m), A(n, m)
      include 'fv.inc'
      integer i, j

      do i = 1, n
      do j = 1, m
        X(i,j) = -A(i,j)
      end do
      end do

      return
      end
      integer function fvFTA(X, F, A, n, m)
C--------------------------------------------------------------------------
C -- calculate X = F*A for n x m matrices X, A and double precision number F
C -- it is allowed for result matrix X to be A
C
      implicit none
      integer n, m
      double precision X(n,m), F, A(n, m)
      include 'fv.inc'
      integer i, j

      fvFTA = NORMAL

      do i = 1, n
      do j = 1, m
        X(i,j) = F*A(i,j)
      end do
      end do

      return
      end
      subroutine fvAPB(X, A, B, n, m)
C--------------------------------------------------------------------------
C -- calculate X = A + B for n x m matrices X, A, B
C
      implicit none
      integer n, m
      double precision X(n,m), A(n, m), B(n, m)
      include 'fv.inc'
      integer i, j

      do i = 1, n
      do j = 1, m
        X(i,j) = A(i,j) + B(i,j)
      end do
      end do

      return
      end
      subroutine fvAPBT(X, A, B, n, m)
C--------------------------------------------------------------------------
C -- calculate X = A + B^T for n x m matrices X, A, B
C -- note: X can be equal to A, but not to B
C
      implicit none
      integer n, m
      double precision X(n,m), A(n, m), B(n, m)
      include 'fv.inc'
      integer i, j

      do i = 1, n
      do j = 1, m
        X(i,j) = A(i,j) + B(j,i)
      end do
      end do

      return
      end
      subroutine fvAB(X, A, B, n0, n1, n2)
C--------------------------------------------------------------------------
C -- calculate X = A.B for n0 x n1 matrix A and n1 x n2 matrix B
C -- 
      implicit none
      integer n0, n1, n2
      double precision X(n0,n2), A(n0,n1), B(n1,n2)
      include 'fv.inc'
      integer i, j, k
      double precision d

      do i = 1, n0
      do k = 1, n2

        d = 0d0
        do j = 1, n1
          d = d + A(i,j)*B(j,k)
        end do
        X(i,k) = d

      end do
      end do

      return
      end
      subroutine fvATB(X, A, B, n0, n1, n2)
C--------------------------------------------------------------------------
C -- calculate n1 x n2 matrix X = A^T.B 
C -- with n0 x n1 matrix A and n0 x n2 matrix B
C -- 
C
      implicit none
      integer n0, n1, n2
      double precision X(n1,n2), A(n0,n1), B(n0,n2)
      include 'fv.inc'
      integer i, k, l
      double precision d

      do i = 1, n1
      do l = 1, n2

        d = 0d0
        do k = 1, n0
          d = d + A(k,i)*B(k,l)
        end do
        X(i,l) = d

      end do
      end do

      return
      end
      subroutine fvABT(X, A, B, n0, n1, n2)
C--------------------------------------------------------------------------
C -- calculate n0 x n2 matrix X = A.B^T
C -- with n0 x n1 matrix A and n2 x n1 matrix B
C -- 
C
      implicit none
      integer n0, n1, n2
      double precision X(n0,n2), A(n0,n1), B(n2,n1)
      include 'fv.inc'
      integer i, k, l
      double precision d

      do i = 1, n1
      do l = 1, n2

        d = 0d0
        do k = 1, n0
          d = d + A(i,k)*B(l,k)
        end do
        X(i,l) = d

      end do
      end do

      return
      end
      subroutine fvATBC(X, A, B, C, n, m)
C--------------------------------------------------------------------------
C -- calculate m x m matrix X = A^T.B.C 
C -- for n x n matrix B
C -- and n x m matrices A,C
C -- 
C
      implicit none
      integer m, n
      double precision X(m,m), A(n,m), B(n,n), C(n,m)
      include 'fv.inc'
      double precision temp(20,20)
      integer i, j, k, l
      double precision d

      do k = 1, n
      do j = 1, m

        d = 0d0
        do l = 1, n
          d = d + B(k,l)*C(l,j)
        end do
        temp(k,j) = d

      end do
      end do

      do i = 1, m
      do j = 1, m

        d = 0d0
        do k = 1, n
          d = d + A(k,i)*temp(k,j)
        end do
        X(i,j) = d

      end do
      end do

      return
      end
      subroutine fvsATBA(X, A,B,n,m)
C--------------------------------------------------------------------------
C -- calculate m x m matrix X = A^T.B.A
C -- for symmetric n x n matrix B
C -- and n x m matrix A
C -- 
C
      implicit none
      integer m, n
      double precision X(m,m), A(n,m), B(n,n)
      include 'fv.inc'
      double precision temp(20,20)
      integer i, j, k, l
      double precision d

      do k = 1, n
      do j = 1, m

        d = 0d0
        do l = 1, n
          d = d + B(k,l)*A(l,j)
        end do
        temp(k,j) = d

      end do
      end do

      do i = 1, m
      do j = 1, i

        d = 0d0
        do k = 1, n
          d = d + A(k,i)*temp(k,j)
        end do
        X(i,j) = d
        if (i.NE.j) X(j,i) = d

      end do
      end do

      do i = 1, m
      do j = 1, i-1

        X(j,i) = X(i,j)

      end do
      end do

      return
      end
      subroutine fvABCT(X, A, B, C, m, n)
C--------------------------------------------------------------------------
C -- calculate m x m matrix X = A.B.C^T 
C -- for n x n matrix B
C -- and m x n matrices A, C
C -- 
C
      implicit none
      integer m, n
      double precision X(m,m), A(m,n), B(n,n), C(m,n)
      include 'fv.inc'
      double precision temp(20,20)
      integer i, j, k, l
      double precision d

      do i = 1, n
      do l = 1, m
        d = 0d0
        do j = 1, n
          d = d + B(i,j)*C(l,j)
        end do
        temp(i,l) = d
      end do
      end do

      do k = 1, m
      do l = 1, m
        d = 0d0
        do i = 1, n
          d = d + A(k,i)*temp(i,l)
        end do
        X(k,l) = d
      end do
      end do



      return
      end
      subroutine fvsABAT(X, A, B, m, n)
C--------------------------------------------------------------------------
C -- calculate m x m matrix X = A.B.A^T 
C -- for a symmetric n x n matrix B and m x n matrix A
C -- 
C
      implicit none
      integer m, n
      double precision X(m,m), A(m,n), B(n,n)
      include 'fv.inc'
      double precision temp(20,20)
      integer i, j, k, l
      double precision d

      do i = 1, n
      do l = 1, m
        d = 0d0
        do j = 1, n
          d = d + B(i,j)*A(l,j)
        end do
        temp(i,l) = d
      end do
      end do

      do k = 1, m
      do l = 1, k
        d = 0d0
        do i = 1, n
          d = d + A(k,i)*temp(i,l)
        end do
        X(k,l) = d
        if (l.NE.k) X(l,k) = d
      end do
      end do

      return
      end
      subroutine fvABTCT(X, A, B, C, m, n)
C--------------------------------------------------------------------------
C -- calculate m x m matrix X = A.B^T.C^T 
C -- for n x n matrix B
C -- and m x n matrices A, C
C -- 
C
      implicit none
      integer m, n
      double precision X(m,m), A(m,n), B(n,n), C(m,n)
      include 'fv.inc'
      double precision temp(20,20)
      integer i, j, k, l
      double precision d

      do i = 1, n
      do l = 1, m
        d = 0d0
        do j = 1, n
          d = d + B(j,i)*C(l,j)
        end do
        temp(i,l) = d
      end do
      end do

      do k = 1, m
      do l = 1, m
        d = 0d0
        do i = 1, n
          d = d + A(k,i)*temp(i,l)
        end do
        X(k,l) = d
      end do
      end do

      return
      end
      double precision function fvProb(chi2, ndof)
C----------------------------------------------------------------------
C -- check for negative chi2 and p
C 
      implicit none
      double precision chi2  
      integer ndof
      include 'fv.inc'
      real PROB

      if (chi2 .le. 0d0) then
        write(text,'(a,g12.4,a)') 
     1      'Error in fvProb: chi**2 = ', chi2, ' <= 0.'
        call fvError(text)
        fvProb = 0d0
      else
        fvProb = PROB(real(chi2), ndof)
      end if

      return
      end
      double precision function fvPyth(a,b)
C--------------------------------------------------------------------
C -- calculate sqrt(a*a + b*b) w/o destructive underflow or overflow
C
      implicit none
      double precision a,b,at,bt,ct
      double precision VERYSMALL
      parameter (VERYSMALL=1d-18)

      at = abs(a)
      bt = abs(b)

      if (at .gt. bt) then
         ct = bt/at
C -- protect against underlow
         if (ct .lt. VERYSMALL) ct = 0d0
         fvPyth = at*sqrt(1d0+ct*ct)
      else
         if (bt .ne. 0d0) then
            ct = at/bt
            if (ct .lt. VERYSMALL) ct = 0d0
            fvPyth = bt*sqrt(1d0+ct*ct)
         else
            fvPyth = 0d0
         end if
      end if

      return
      end
      integer function fvSVD(a,m,n,mp,np,w,v)
C--------------------------------------------------------------------
C -- do a singular value decomposition of matrix A
C -- 
C -- origin: W.H. Press et al., Numerical Recipes, Cambridge Univ. Press
C -- 
C -- Given a matrix A, with logical dimensions m by n
C -- and physical dimensions mp by np, this
C -- routine computes its singular value decomposition 
C -- A=U.W.V^T. The matrix U replaces
C -- A on output. The diagonal matrix of singular values W
C -- is output as a vector W. The matrix
C -- V (not the transpose V^T) is output as V.
C -- m must be greater or equal to n; if it is smaller
C -- then A should be filled up to square with zero rows.
C--
      implicit none

      integer m,n,mp,np
      double precision a(mp,np), w(np),v(np,np)

C -- maximum anticipated value of n
      integer nmax
      parameter (nmax=100)
      double precision rv1(nmax)
      double precision g,scale,anorm,s,f,h,x,y,z,c
      integer i,j,k,l,jj,nm,its

      double precision fvPyth
      external fvPyth

      integer NORMAL, ERROR
      parameter (NORMAL = 1, ERROR = 2)


      fvSVD = NORMAL

      if (m.lt.n) then
         fvSVD = ERROR
         return
      end if

C -- Householder reduction to bidiagonal form
      g=0d0
      scale=0d0
      anorm=0d0
      do i=1,n
         l=i+1
         rv1(i)=scale*g
         g=0d0
         s=0d0
         scale=0d0
         if (i.le.m) then
            do k=i,m
               scale=scale+abs(a(k,i))
            end do
            if (scale .ne. 0d0) then
               do k=i,m
                  a(k,i)=a(k,i)/scale
                  s=s+a(k,i)*a(k,i)
               end do
               f=a(i,i)
               g=-sign(sqrt(s),f)
               h=f*g-s
               a(i,i)=f-g
               if (i.ne.n) then
                  do j=l,n
                     s=0d0
                     do k=i,m
                        s=s+a(k,i)*a(k,j)
                     end do
                     f=s/h
                     do k=i,m
                        a(k,j)=a(k,j)+f*a(k,i)
                     end do
                  end do
               end if
               do k=i,m
                  a(k,i)=scale*a(k,i)
               end do
            end if
         end if
         w(i)=scale*g
         g=0d0
         s=0d0
         scale=0d0
         if ((i.le.m).and.(i.ne.n)) then
            do k=l,n
               scale=scale+abs(a(i,k))
            end do
            if (scale .ne. 0d0) then
               do k=l,n
                  a(i,k)=a(i,k)/scale
                  s=s+a(i,k)*a(i,k)
               end do
               f=a(i,l)
               g=-sign(sqrt(s),f)
               h=f*g-s
               a(i,l)=f-g
               do k=l,n
                  rv1(k)=a(i,k)/h
               end do
               if (i.ne.m) then
                  do j=l,m
                     s=0d0
                     do k=l,n
                        s=s+a(j,k)*a(i,k)
                     end do
                     do k=l,n
                        a(j,k)=a(j,k)+s*rv1(k)
                     end do
                  end do
               end if
               do k=l,n
                  a(i,k)=scale*a(i,k)
               end do
            end if
         end if
         anorm=max(anorm,(abs(w(i))+abs(rv1(i))))
      end do
C -- Accumulation of right-hand transformation
      do i=n,1,-1
         if (i.lt.n) then
            if (g.ne.0d0) then
               do j=l,n
C -- double division to avoid possible underflow
                  v(j,i)=(a(i,j)/a(i,l))/g
               end do
               do j=l,n
                  s=0d0
                  do k=l,n
                     s=s+a(i,k)*v(k,j)
                  end do
                  do k=l,n
                     v(k,j)=v(k,j)+s*v(k,i)
                  end do
               end do
            end if
            do j=l,n
               v(i,j)=0d0
               v(j,i)=0d0
            end do
         end if
         v(i,i)=1d0
         g=rv1(i)
         l=i
      end do
C -- Accumulation of left-hand transformation
      do i=n,1,-1
         l=i+1
         g=w(i)
         if (i.lt.n) then
            do j=l,n
               a(i,j)=0d0
            end do
         end if
         if (g.ne.0d0) then
            g=1d0/g
            if (i.ne.n) then
               do j=l,n
                  s=0d0
                  do k=l,m
                     s=s+a(k,i)*a(k,j)
                  end do
                  f=(s/a(i,i))*g
                  do k=i,m
                     a(k,j)=a(k,j)+f*a(k,i)
                  end do
               end do
            end if
            do j=i,m
               a(j,i)=a(j,i)*g
            end do
         else
            do j=i,m
               a(j,i)=0d0
            end do
         end if
         a(i,i)=a(i,i)+1d0
      end do
C -- Diagonalization of the bidiagonal form
      do k=n,1,-1
C -- loop over singular values
         do its=1,30
C -- loop over allowed iterations
            do l=k,1,-1
C -- Test for splitting
               nm=l-1
C -- note tha rv1(1) is always zero
               if ((abs(rv1(l))+anorm).eq.anorm) goto 2
               if ((abs(w(nm))+anorm).eq.anorm) goto 1
            end do
 1          continue
            c=0d0
C -- Cancellation of rv1(l) if l>1
            s=1d0
            do i=l,k
               f=s*rv1(i)
               if ((abs(f)+anorm).ne.anorm) then
                  g=w(i)
                  h=fvPyth(f,g)
                  w(i)=h
                  h=1d0/h
                  c= (g*h)
                  s=-(f*h)
                  do j=1,m
                     y=a(j,nm)
                     z=a(j,i)
                     a(j,nm)=(y*c)+(z*s)
                     a(j,i)=-(y*s)+(z*c)
                  end do
               end if
            end do
 2          continue
            z=w(k)
            if (l.eq.k) then
C -- Convergence
               if (z.lt.0d0) then
C -- singular value is made nonnegative
                  w(k)=-z
                  do j=1,n
                     v(j,k)=-v(j,k)
                  end do
               end if
               goto 3
            end if
            if (its.eq.30) then
CCCCC               call fvError('No convergence in 30 iterations')
               fvSVD = ERROR
               return
            end if
            x=w(l)
C -- shift from bottom 2-by-2 mirror
            nm=k-1
            y=w(nm)
            g=rv1(nm)
            h=rv1(k)
            f=((y-z)*(y+z)+(g-h)*(g+h))/(2d0*h*y)
            g=fvPyth(f,1d0)
            f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
C -- next QR transformation
            c=1d0
            s=1d0
            do j=l,nm
               i=j+1
               g=rv1(i)
               y=w(i)
               h=s*g
               g=c*g
               z=fvPyth(f,h)
               rv1(j)=z
               c=f/z
               s=h/z
               f= (x*c)+(g*s)
               g=-(x*s)+(g*c)
               h=y*s
               y=y*c
               do jj=1,n
                  x=v(jj,j)
                  z=v(jj,i)
                  v(jj,j)= (x*c)+(z*s)
                  v(jj,i)=-(x*s)+(z*c)
               end do
               z=fvPyth(f,h)
               w(j)=z
C -- rotation can be arbitrary if z=0
               if (z.ne.0d0) then
                  z=1d0/z
                  c=f*z
                  s=h*z
               end if
               f= (c*g)+(s*y)
               x=-(s*g)+(c*y)
               do jj=1,m
                  y=a(jj,j)
                  z=a(jj,i)
                  a(jj,j)= (y*c)+(z*s)
                  a(jj,i)=-(y*s)+(z*c)
               end do
            end do
            rv1(l)=0d0
            rv1(k)=f
            w(k)=x
         end do
 3       continue
      end do
      return
      end
      integer function fvLUdcmp(a,n,np,indx,d)
C----------------------------------------------------------------------- 
C -- do a LU-decomposition
C -- origin: W.H. Press et al., Numerical Recipes, Cambridge Univ. Press
C --
C -- given an NxN matrix A, with physical dimensions NP, 
C -- this routine replaces it by the LU
C -- decomposition of a rowwise permutation of itself.
C -- A and N are input. A is output, arranged
C -- as in equation (2.3.14) above; INDX is an output
C -- vector which records the row permutation
C -- effected by the partial pivoting; D is output
C -- as +- 1 depending on whatever the numer of
C -- row interchanges was even or odd, respectively.
C -- This routine is used in combination with
C -- LUBKSB to solve linear equations or invert a matrix.
C
      implicit none
      integer n,np
      double precision a(np,np),indx(n),d

      integer ERROR,NORMAL
      parameter (ERROR=2,NORMAL=1)
C -- largest expected N, a small number
      integer NMAX
      double precision TINY
      parameter (NMAX=100, TINY = 1.0d-30)

C -- VV stores the implicit scaling of each row
      double precision vv(NMAX),aamax,sum,dum
      integer i,j,k,imax

      fvLUdcmp = NORMAL

      d=1d0
      do i=1,n
         aamax=0d0
         do j=1,n
            if (abs(a(i,j)) .GT. aamax) aamax=abs(a(i,j))
         end do
         if (aamax .EQ. 0d0) then
C -- matrix singular, no nonzero largest element
            fvLUdcmp=ERROR
            return
         end if
C -- save the scaling
         vv(i)=1d0/aamax
      end do
      do j=1,n
C -- loop over columns of Crout's method
         do i=1,j-1
C -- this is equation (2.3.12) except for i=j
            sum=a(i,j)
            do k=1,i-1
               sum=sum-a(i,k)*a(k,j)
            end do
            a(i,j)=sum
         end do
C -- initialize for the search for the largest pivot element
         aamax=0d0
         do i=j,n
C -- this is i=j of equation (2.3.12) and i=j+1..n of equation (2.3.13)
            sum=a(i,j)
            do k=1,j-1
               sum=sum-a(i,k)*a(k,j)
            end do
            a(i,j)=sum
C -- figure of merit for the pivot
            dum=vv(i)*abs(sum)
C -- is it better than the best so far?
            if (dum .GE. aamax) then
               imax=i
               aamax=dum
            endif
         end do
C -- do we need to interchange rows?
         if (j .NE. imax) then
C -- yes, do so
            do k=1,n
               dum=a(imax,k)
               a(imax,k)=a(j,k)
               a(j,k)=dum
            end do
C -- change the parity of d
            d=-d
C -- also interchange the scale factor
            vv(imax)=vv(j)
         endif
         indx(j)=imax
C -- now, finally, divide by the pivot element
C -- if the pivot element is zero the matrix is singulax
C -- (al leas to the precsision of the algorithm)
C -- for some applications on singular matrices it is desirable to 
C -- substitute TINY for zero
         if (a(j,j) .EQ. 0d0) a(j,j) = TINY
         if (j .NE. N) then
            dum=1d0/a(j,j)
            do i=j+1,n
               a(i,j)=a(i,j)*dum
            end do
         end if
C -- go back to the next column in the reduction
      end do

      return
      end
      subroutine fvLUbksb(a,n,np,indx,b)
C----------------------------------------------------------------------- 
C -- do a LU back substitution
C -- origin: W.H. Press et al., Numerical Recipes, Cambridge Univ. Press
C --
C -- Solves the set of N linear equations A.X=B. 
C -- Here A is input, not as the matrix A
C -- but rather as its LU decomposition,
C -- determined by the routine fvLUdcmp. INDX is input
C -- as the permutation vector returned by fvLUdcmp.
C -- B is input as the right-hand side vector
C -- B, and returns with the solution vector X.
C -- A, N, NP, and INDX are not modified by this
C -- routine and can be left in place for successive calls
C -- with different right-hand sides B. This
C -- routine takes into account the possibility that 
C -- B will begin with many zero elements, so
C -- it is effiscient for use in matrix inversion.
C
      implicit none
      integer n,np
      double precision a(np,np),indx(n),b(n)
      integer ii,i,j,k,ll
      double precision sum

C -- when ii is set to a positive value, it will become the index of the first
C -- nonvanishing element of b. We now do the forward substitution,
C -- equation (2.3.6). The only new wrinkle is to unscramble
C -- the permutation as we go
      ii=0
      do i=1,n
         ll=indx(i)
         sum=b(ll)
         b(ll)=b(i)
         if (ii .NE. 0) then
            do j=ii,i-1
               sum=sum-a(i,j)*b(j)
            end do
         else if (sum .NE. 0d0) then
C -- a nonzero element was encountered, so from now on we will have to
C -- do the sums in the loop above
            ii=i
         end if
         b(i)=sum
      end do
C -- now we do the backsubstitution, equation 2.3.7
      do i=n,1,-1
         sum=b(i)
         if (i.LT.n) then
            do j=i+1,n
               sum=sum-a(i,j)*b(j)
            end do
         end if
C -- store a component of the solution vector X
         b(i)=sum/a(i,i)
      end do
C -- all done!
      return
      end
      subroutine fvTinv(a,b,n)
C-------------------------------------------------------------------------
C -- test inversion of a
C -- a = b^(-1)
      implicit none
      integer n
      double precision a(n,n),b(n,n)

      integer i0,i1,j
      integer NMAX
      parameter (NMAX=6)
      double precision c(NMAX,NMAX), sum
      

      do i0 = 1, n
         do i1 =1, n
            sum = 0d0
            do j=1,n
               sum = sum + a(i0,j)*b(j,i1)
            end do
            c(i0,i1) = sum
         end do
         if (a(i0,i0) .LT. 0d0) then
            print *,'ERROR: error matrix a has neg. entries on diag'
         end if
         if (b(i0,i0) .LT. 0d0) then
            print *,'ERROR: error matrix b has neg. entries on diag'
         end if
      end do
      do i0=1,n
         sum = 0d0
         do i1=1,n
            if (i0 .ne. i1) sum = sum+a(i0,i1)
         end do
         print *, i0, a(i0,i0), sum/dble(n)
      end do
      return
      end
