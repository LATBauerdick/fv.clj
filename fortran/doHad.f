      subroutine doHad
C --------------------------------------------------------------------
      implicit none
      integer KRUN, KEVT
      integer NEVT
      parameter (NEVT=1)
      integer nTuple
      parameter (nTuple = 6)
      integer nTuple2
      parameter (nTuple2 = 10)
      integer PAWLEN
      parameter (PAWLEN = 100000)


      include 'fvt.inc'
      integer runList(2,NEVT)
      data runList /
     1  9096,2
     1/
      integer trList(40,0:0)
      data trList /
     1  1,2,3,4,5,6,7,8,9,10,
     1  11,12,13,14,15,16,17,18,19,20,
     1  21,22,23,24,25,26,27,28,29,30,
     1  31,32,33,34,35,36,37,38,39,40    
     1/
      integer pList(40,0:0)
C -- 1:electron, 2:muon, 3:tau, 4:pi
      data pList /
     1  4,4,4,4,4,4,4,4,4,4,
     1  4,4,4,4,4,4,4,4,4,4,
     1  4,4,4,4,4,4,4,4,4,4,
     1  4,4,4,4,4,4,4,4,4,4    
     1/
      integer i, j, k, it, i0, i1
      integer status

      double precision h(5,MAXTR), ChLam(5,MAXTR), ChU(5,5,MAXTR), x(3)
      double precision chi2, p(3),Dp(3,3),Ep(3,3),p0(3),pp0(3),
     1                 v00(3),v0(3),Cv0(3,3),tau,chi2tau,pp,sigTau
      logical first
      data first/.true./
      SAVE first
      integer nMC
      data nMC/1000/

      real PROB, RG32
      external PROB, RG32
      
      real nt(nTuple)
      character*8 nTags(nTuple)
      data nTags /
     1  'prob',
     1  'pw',
     1  'ptl',
     1  'ppsi0',
     1  'pd0',
     1  'pz0'
     1/
      real nt2(nTuple2)
      character*8 nTags2(nTuple2)
      data nTags2 /
     1  'prob',
     1  'vx',
     1  'vy',
     1  'vz',
     1  'px',
     1  'py',
     1  'pz',
     1  'tau',
     1  'sigTau',
     1  'probTau'
     1/
      integer HMEMOR(PAWLEN)
      common /PAWC/ HMEMOR

      integer fvtRead, fvtDo
      external fvtRead, fvtDo

      external HBOOKN, HROPEN
      integer LRECL, ISTAT,ICYCLE

      call HLIMIT(PAWLEN)
      LRECL=1024
      call HROPEN(48, 'NT', 'had-0.nt', 'N', LRECL, ISTAT)
      call HBOOKN(1, 'tracks', nTuple, 'NT', 1000, nTags)
      call HBOOKN(2, 'vertex', nTuple2, 'NT', 1000, nTags2)



      do i = 1, NEVT
        KRUN = runList(1,i)
        KEVT = runList(2,i)

        if (fvtPrint) then
          write(fvtPlun,'(1x,78(1H=))') 
          write(fvtPlun,'(/1x,a,i5,a,i5)') 
     1      'Run ', KRUN, ' Evt ', KEVT
        end if
C -- read event from private data file
        status = fvtRead(KRUN,KEVT)
C -- set vertex to approx zero
        do i1 = 1,3
          tx(i1,1) = 0.0001d0
        end do
C -- get private copy of vertex
        call fvCopy(x, tx(1,1),3,1)
C -- fill error matrix acc. to beam spot size
        do i0 = 1,3
        do i1 = 1,3
          tCx(i0,i1,1) = 0.0d0
        end do
        end do
        tCx(1,1,1) = 0.04d0**2   
        tCx(2,2,1) = 0.004d0**2   
        tCx(3,3,1) = 4.d0**2      

C -- diagonalize covariance matrices of helix parameters
        do j = 1, tnt
          call fvtMCdiag(ChLam(1,j),ChU(1,1,j), tCh(1,1,j))
        end do
C -- force all tracks to come from common vertex
        do j = 1, tnt
          call fvtMCunify(h(1,j), th(1,j))
        end do

C -- set point v0 and covariance matrix
        do i0 = 1,3
        do i1 = 1,3
          Cv0(i0,i1) = 0.0d0
        end do
        end do
        Cv0(1,1) = 0.04d0**2   
        Cv0(2,2) = 0.004d0**2  
        Cv0(3,3) = 4.d0**2   
        tau = 1.
        
        do j = 1, nMC
C -- smear inital vertex according to tCx(i,j,1)
C -- (assmued to be diagonal)
          do k = 1, 3
            tx(k,1) = x(k) + sqrt(tCx(k,k,1))*RG32(k)
          end do
CCC -- if no smear:
CC            call fvCopy(tx(1,1),x,3,1)
C -- smear track parameters according to cov. matrix
          do k = 1, tnt
            it = trList(k,0)
            call fvtMCsmear(th(1,it), h(1,it),ChLam(1,it),ChU(1,1,it))
CCC -- if no smear:
CC            call fvCopy(th(1,it), h(1,it),5,1)
          end do

C -- do the vertex fit
          status = fvtDo(chi2, tnt,trList(1,0))

C -- sum up momentum of half of the tracks
          status = fvSumQ(p,Dp,Ep,
     1                tnt+1,(tnt+1)/2,trList(1,0),tw2pt,
     1                tx(1,1),tCx(1,1,1),tq,tCq,th,tGh)
C -- fit distance to point v0, tau = 1.5d0
          if (first) then
            pp = sqrt(p(1)**2+p(2)**2+p(3)**2)
            v0(1) = x(1) - p(1)/pp*1.5d0
            v0(2) = x(2) - p(2)/pp*1.5d0
            v0(3) = x(3) - p(3)/pp*1.5d0
            first = .false.
          end if
          do k = 1, 3
            v00(k) = v0(k)+Cv0(k,k)*RG32(k)
          end do
          status = fvFitD(tau,sigTau,chi2tau,
     1                tau,v00,Cv0,tx(1,1),tCx(1,1,1),p,Dp,Ep)

C -- calculate and sum up unsmeared momenta from unsmeared helix param's
          call fvZeroA(p0,3,1)
          do k = 1, (tnt+1)/2
            it = trList(k,0)
            call fvq2pvec(pp0, h(1,it),tw2pt)
            p0(1) = p0(1)+pp0(1)
            p0(2) = p0(2)+pp0(2)
            p0(3) = p0(3)+pp0(3)
          end do

C -- histogram the results
          do k = 1, tnt
            it = trList(k,0)
            nt(1) = tChi2(it)
            nt(1) = PROB(nt(1),2)

            nt(2) = (th(1,it)-h(1,it))/sqrt(tCh(1,1,it))
            nt(3) = (th(2,it)-h(2,it))/sqrt(tCh(2,2,it))
            nt(4) = (th(3,it)-h(3,it))/sqrt(tCh(3,3,it))
            nt(5) = (th(4,it)-h(4,it))/sqrt(tCh(4,4,it))
            nt(6) = (th(5,it)-h(5,it))/sqrt(tCh(5,5,it))
            call HFN(1, nt)
          end do

          nt2(1) = chi2
          nt2(1) = PROB(nt2(1),2*tnt)
          nt2(2) = tx(1,1)/sqrt(tCx(1,1,1))
          nt2(3) = tx(2,1)/sqrt(tCx(2,2,1))
          nt2(4) = tx(3,1)/sqrt(tCx(3,3,1))
          nt2(5) = (p(1)-p0(1))/sqrt(Dp(1,1))
          nt2(6) = (p(2)-p0(2))/sqrt(Dp(2,2))
          nt2(7) = (p(3)-p0(3))/sqrt(Dp(3,3))
          nt2(8) = tau
          nt2(9) = sigTau
          nt2(10) = chi2tau
          nt2(10) = PROB(nt2(10),1)
          call HFN(2, nt2)
        end do

      end do

      call HROUT(1, ICYCLE, ' ')
      call HROUT(2, ICYCLE, ' ')
      call HREND('NT')

      return
      end
