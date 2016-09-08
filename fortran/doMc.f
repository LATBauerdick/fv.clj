      subroutine doMc
C --------------------------------------------------------------------
      implicit none
      include 'fvt.inc'

      integer NEVT, NTRK
      parameter (NEVT=1,
     1           NTRK=4)
      integer runList(2,NEVT)
      data runList / 
     1  5166,000981
     1/
C -- tracks:
C -- event: Z0 -> l0 l2
C --              l0 -> l1 V+ V-
C -- V+, V-, l1, l2
      integer tList(NTRK,0:NEVT)
      data tList /
     1  1,2,3,4,    
     1  4,3,2,1
     1/
      integer pList(NTRK,0:NEVT)
C -- 1:electron, 2:muon, 3:tau, 4:pi
      data pList /
     1  4,4,4,4,    
     1  2,2,4,2
     1/
      double precision mass(4)
      data mass/
     1  0.00051099906d0,
     1  0.105658387d0,
     1  1.7841d0,
     1  0.1395675d0
     1/

C -- HBOOK related...
      integer PAWLEN
      parameter (PAWLEN = 100000)
      integer nTuple
      parameter (nTuple = 22)
      real nt1(nTuple)
      character*8 nTags(nTuple)
      data nTags /
     1  'undef',
     1  'd0',
     1  'sigd0',
     1  'd',
     1  'sigd',
     1  'prd',
     1  'm',
     1  'sigm',
     1  'vpx',    
     1  'vpy',    
     1  'vpz',    
     1  'sigvpx',
     1  'sigvpy',
     1  'sigvpz',
     1  'prp', 
     1  'vsx',    
     1  'vsy',    
     1  'vsz',    
     1  'sigvsx',
     1  'sigvsy',
     1  'sigvsz',
     1  'prs'    
     1/
      integer HMEMOR(PAWLEN)
      common /PAWC/ HMEMOR
      external HBOOKN, HROPEN
      integer LRECL, ISTAT,ICYCLE
C -- end of HBOOK stuff

      integer i, j, it, new, status
      integer KRUN, KEVT

      double precision h(5,MAXTR), ChLam(5,MAXTR), ChU(5,5,MAXTR)
      
      double precision p(3),Dp(3,3),Ep(3,3),v(3),Cv(3,3)
      double precision vp(3),vs(3)
      double precision chi2, pr
      double precision tau, tauNom, m, sigM, ml(MAXTR),
     1                 d0,sigd0,d,sigd,chi2d,prd,
     1                 chi2p,prp
      integer tpList(2)
      data tauNom /1d0/

      double precision ps(3),psAbs, p0(3),p1(3),q(3)

      double precision bs(3), Cbs(3,3)
      data Cbs/0.0016d0,0d0,0d0,
     1     0d0,1.6d-5,0d0,
     1     0d0,0d0,4
     1/
      integer imc, nmc
      data nmc/1000/
      
      real RN32, RG32
      external RN32, RG32
      integer fvtRead, fvtDo
      external fvtRead, fvtDo

C -- init HBOOK
      call HLIMIT(PAWLEN)
      LRECL=1024
      call HROPEN(48, 'NT', 'mc-0.nt', 'N', LRECL, ISTAT)
      call HBOOKN(1, 'mc', nTuple, 'NT', 1000, nTags)

C -- init some variables

C -- beam spot location
      bs(1) = sqrt(Cbs(1,1))*RG32(bs(1))
      bs(2) = sqrt(Cbs(2,2))*RG32(bs(2))
      bs(3) = sqrt(Cbs(3,3))*RG32(bs(3))

      do i = 1, NEVT
        KRUN = runList(1,i)
        KEVT = runList(2,i)
        write(fvtPlun,'(1x,78(1H=))') 
        write(fvtPlun,'(/1x,a,i5,a,i5)') 
     1    'Run ', KRUN, ' Evt ', KEVT

C -- get helix parameters of event
        status = fvtRead(KRUN,KEVT)

C
C --
C -- setup MC truth input based on real event topology:
C --
C -- generate a primary vertex inside the beamspot
        vp(1) = bs(1)+sqrt(Cbs(1,1))*RG32(vp(1))
        vp(2) = bs(2)+sqrt(Cbs(2,2))*RG32(vp(2))
        vp(3) = bs(3)+sqrt(Cbs(3,3))*RG32(vp(3))
C -- move l+l- to primary vertex 
        call fvZeroA(v,3,1)
        do j = 3, 4
          it = tList(j,i)
          call fvq(q, v,th(1,it))
          call fvh(th(1,it), vp,q)
        end do

C -- calculate secondary vertex in direction of flight of V+V-
        call fvq(q, v,th(1,tList(1,i)))
        call fvp3(p0, q,tw2pt)

        call fvq(q, v,th(1,tList(2,i)))
        call fvp3(p1, q,tw2pt)

        ps(1) = p0(1)+p1(1)
        ps(2) = p0(2)+p1(2)
        ps(3) = p0(3)+p1(3)
        psAbs = sqrt(ps(1)**2+ps(2)**2+ps(3)**2)

 1      tau = RN32(tau)
        if (tau .EQ. 0.) goto 1
Cdon't smear tau at the moment        tau = -log(tau/tauNom)
        tau = tauNom

        vs(1) = vp(1)+tau*ps(1)/psAbs
        vs(2) = vp(2)+tau*ps(2)/psAbs
        vs(3) = vp(3)+tau*ps(3)/psAbs

C -- move V+V- pair to secondary vertex vs
        do j = 1, 2
          it = tList(j,i)
          call fvq(q, v,th(1,it))
          call fvh(th(1,it), vs,q)
        end do

C -- store a copy of un-smeared tracks
C -- diagonalize helix covariance matrices for all tracks
C -- for MC smearing of helix parameters
        do j = 1, tnt
          call fvCopy(h(1,j), th(1,j), 5,1)
          call fvtMCdiag(ChLam(1,j),ChU(1,1,j), tCh(1,1,j))
        end do

C -- do MC 
        do imc = 1, nmc
C -- smear helix parameters according to error matrices
C -- using diagonalization parameters ChLam and ChU
          do j = 1, tnt
            call fvtMCsmear(th(1,j), h(1,j),ChLam(1,j),ChU(1,1,j))
CCC -- if no smear: copy from private copy (original)
CC            call fvCopy(th(1,j), h(1,j),5,1)
          end do

C -- do vertex fit to V+V- and store fitted vertex in v, Cv
C -- large initial error matrix for first guess matrix
          call fvFTA(tCx(1,1,1), 10 000d0, tCx(1,1,1),3,3)
          status = fvtDo(chi2, 2,tList(1,i))
          pr = fvProb(chi2,2*2-3)
          call fvCopy(v, tx,3,1)
          call fvCopy(Cv, tCx,3,3)

C -- calculate invariant mass of V+V-
          ml(tList(1,i)) = mass(pList(1,i))
          ml(tList(2,i)) = mass(pList(2,i))
          call fvInvM(m,sigM,
     1         2,tList(1,i),tw2pt,ml,v,Cv,tq,tCq,tGh)

C -- sum up momenta of V+V- and put resultant track parameters into th/tGh
          new = tnt+1
          status = fvSumQ(p,Dp,Ep,
     1                new,2,tList(1,i),tw2pt,v,Cv,tq,tCq,th,tGh)

          d0 = tauNom
C -- calculate distance between vertex of V+V- and beam spot
          call fvDist(d0,sigd0, 
     1                bs,Cbs,v,Cv)

C -- fit distance between the two
          status = fvFitD(d,sigd,chi2d,
     1                    d,bs,Cbs,v,Cv,p,Dp,Ep)
          if (iand(status,1).NE.1) then
            print *, 'fit to distance did not work...'
          end if
          prd = fvProb(chi2d, 1)

C -- fit vertex of sum particle and l(V)
          tpList(1) = tList(3,i)
          tpList(2) = new
C -- large initial error matrix for first guess matrix
          call fvFTA(tCx(1,1,1), 10 000d0, tCx(1,1,1),3,3)
          status = fvtDo(chi2p, 2,tpList(1))
          prp = fvProb(chi2p,2*2-3)
          
C -- fill ntuple
          nt1(2) = d0
          nt1(3) = sigd0
          nt1(4) = d
          nt1(5) = sigd
          nt1(6) = prd
          nt1(7) = m
          nt1(8) = sigm
C -- reconstructed primary vertex
          nt1(9) = tx(1,1) - vp(1)
          nt1(10) = tx(2,1) - vp(2)
          nt1(11) = tx(3,1) - vp(3)
          nt1(12) = sqrt(tCx(1,1,1))
          nt1(13) = sqrt(tCx(2,2,1))
          nt1(14) = sqrt(tCx(3,3,1))
          nt1(15) = prp

C -- reconstructed secondary  vertex
          nt1(16) = v(1) - vs(1)
          nt1(17) = v(2) - vs(2)
          nt1(18) = v(3) - vs(3)
          nt1(19) = sqrt(Cv(1,1))
          nt1(20) = sqrt(Cv(2,2))
          nt1(21) = sqrt(Cv(3,3))
          nt1(22) = pr

          call HFN(1, nt1)

        end do
      end do

      call HROUT(1, ICYCLE, ' ')
      call HREND('NT')

      return
      end
