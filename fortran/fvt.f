C----------------------------------------------------------------------
C
C FVT, test the FV package to fit vertices and find secondary vertices
C
C author: Lothar A.T. Bauerdick, CERN/PPE
C date:   Jan. 1991
C
C ref.:   R. Fruehwirt, ``Applications of Filter Methods...''
C         HEPHY-PUB 516/88, Vienna
C
C----------------------------------------------------------------------
      program fvt
      implicit none
      include 'fvt.inc'
      character*132 text

      integer m5, m5mc, had, wu, mc
      parameter (m5=1, m5mc=2, had=3, wu=4, mc=5)
      integer do
      
      do = m5
      fvtPrint = .true.
      fvtPlun = 6

      call fvSPrt(fvtPrint,fvtPlun)

C -- default value for mag.field
      tw2pt = 4.4968501d-03
      
      call errInit()
      
      if (fvtPrint) then
CLATB just write to stdout      	open(fvtPlun, FILE='fvt.output', STATUS = 'UNKNOWN')
      end if

      if (do .eq. mc) then
CCCCCCCCCCCCCC        call doMc
      else if (do .eq. m5) then
        call doM5
      else if (do .eq. m5mc) then
        call doM5mc
      else if (do .eq. had) then
CCCCCCCCCCCCCC        call doHad
      else if (do .eq. wu) then
        call doWu
      else
        call doOrs
      end if

      call errSumm()

      stop
      end
      subroutine doM5mc
C --------------------------------------------------------------------
      implicit none
      integer KRUN, KEVT
      integer NEVT
      parameter (NEVT=3)
      include 'fvt.inc'
      integer i, i0, i1, status
      double precision m,sigm
      double precision chi2
      integer runList(2,NEVT)
      integer trList(6,0:NEVT)
      integer pList(6,0:NEVT)
      integer fvtRead, fvtDo, fvtInvM
      external fvtRead, fvtDo, fvtInvM

      data runList /
     1  101,7076,
     1  101,8340,
     1  101,12040
     1/
      data trList /
     1  1,2,3,4,5,6,    
     1  1,3,4,5,6,2,
     1  1,2,3,5,6,4,
     1  1,3,4,5,6,2
     1/
C -- 1:electron, 2:muon, 3:tau, 4:pi
      data pList /
     1  4,4,4,4,4,4,    
     1  4,4,4,4,4,4,
     1  4,4,4,4,4,4,
     1  4,4,4,4,4,4
     1/

      do i = 1, NEVT
        KRUN = runList(1,i)
        KEVT = runList(2,i)


        if (fvtPrint) then
          write(fvtPlun,'(1x,78(1H=))') 
          write(fvtPlun,'(/1x,a,i5,a,i5)') 
     1      'Run ', KRUN, ' Evt ', KEVT
        end if
        status = fvtRead(KRUN,KEVT)
        do i0 = 1,3
        do i1 = 1,3
          tCx(i0,i1,1) = tCx(i0,i1,1) * 10 000d0
        end do
        end do
        status = fvtDo(chi2, 6,trList(1,i))

C -- 
        status = fvtRead(KRUN,KEVT)
        do i0 = 1,3
        do i1 = 1,3
          tCx(i0,i1,1) = tCx(i0,i1,1) * 10 000d0
        end do
        end do

        if (fvtPrint) then
          write(fvtPlun,'(/1x,a,i5,a,i5)') 
     1      'common vertex for 5 tracks for Run ', KRUN, 
     1      ' Evt ', KEVT

C -- fill tt, tCt from th, tCh for particles trList, pList
C -- and calculate inv. mass from  non-fitted momenta
          call fvtHelix2P4(tt,tCt,
     1      5,trList(1,i),pList(1,i),th,tCh)
          status = fvtInvM(m,sigm,5,trList(1,i),tt,tCt)
          write(fvtPlun,'(1x,a,g12.5,a,g12.5)') 
     1      'invariant mass is ', m*1000, ' +/- ', sigm*1000
        end if

        status = fvtDo(chi2, 5,trList(1,i))

        if (fvtPrint) then
C -- calculate Energy for particle types in pList
          call fvtQM2P4(tt,tCt, 5,trList(1,i),pList(1,i),tq,tCq)
          status = fvtInvM(m,sigm,5,trList(1,i),tt,tCt)
          write(fvtPlun,'(1x,a,g12.5,a,g12.5)') 
     1      'invariant mass is ', m*1000, ' +/- ', sigm*1000
        end if

      end do

      return
      end
      subroutine doM5
C --------------------------------------------------------------------
      implicit none
      integer KRUN, KEVT
      integer NEVT
      parameter (NEVT=8)
      include 'fvt.inc'
      integer i, i0, i1, status
      double precision m,sigm
      double precision chi2
      integer fvtRead, fvtDo, fvtInvM
      external fvtRead, fvtDo, fvtInvM

      integer runList(2,NEVT)
      integer trList(6,0:NEVT)
      integer pList(6,0:NEVT)
      data runList /
     1  5129,1412,
     1  5158,4656,
     1  5166,1984,
     1  5343,2291,
     1  7849,7984,
     1  8489,4451,
     1  8985,5486,
     1  9019,4769
     1/
      data trList /
     1  1,2,3,4,5,6,    
     1  1,3,4,5,6,2,
     1  1,2,3,5,6,4,
     1  2,3,4,5,6,1,
     1  1,3,4,5,6,2,
     1  1,3,4,5,6,2,
     1  1,2,3,5,6,4,
     1  1,2,3,4,6,5,
     1  1,2,3,4,6,5
     1/
C -- 1:electron, 2:muon, 3:tau, 4:pi
      data pList /
     1  4,4,4,4,4,4,    
     1  4,4,4,4,4,4,
     1  4,4,4,4,4,4,
     1  4,4,4,4,4,4,
     1  4,4,4,4,4,4,
     1  4,4,4,4,4,4,
     1  4,4,4,4,4,4,
     1  4,4,4,4,4,4,
     1  4,4,4,4,4,4
     1/
      do i = 1, NEVT
        KRUN = runList(1,i)
        KEVT = runList(2,i)

        if (fvtPrint) then
          write(fvtPlun,'(1x,78(1H=))') 
          write(fvtPlun,'(/1x,a,i5,a,i5)') 
     1      'Run ', KRUN, ' Evt ', KEVT
        end if

        status = fvtRead(KRUN,KEVT)
        do i0 = 1,3
        do i1 = 1,3
          tCx(i0,i1,1) = tCx(i0,i1,1) * 10 000d0
        end do
        end do
        status = fvtDo(chi2, 6,trList(1,i))

        status = fvtRead(KRUN,KEVT)
        do i0 = 1,3
        do i1 = 1,3
          tCx(i0,i1,1) = tCx(i0,i1,1) * 10 000d0
        end do
        end do

        if (fvtPrint) then
          write(fvtPlun,'(/1x,a,i5,a,i5)') 
     1      'common vertex for 5 tracks for Run ', KRUN, 
     1      ' Evt ', KEVT
C -- fill tt, tCt from th, tCh for particles trList, pList
          call fvtHelix2P4(tt,tCt,
     1      5,trList(1,i),pList(1,i),th,tCh)
C -- calculate inv. mass
          status = fvtInvM(m,sigm,5,trList(1,i),tt,tCt)
          write(fvtPlun,'(1x,a,g12.5,a,g12.5)') 
     1      'invariant mass is ', m*1000, ' +/- ', sigm*1000
        end if

        status = fvtDo(chi2, 5,trList(1,i))

        if (fvtPrint) then
C -- calculate Energy for particle types in pList
          call fvtQM2P4(tt,tCt, 5,trList(1,i),pList(1,i),tq,tCq)
          status = fvtInvM(m,sigm,5,trList(1,i),tt,tCt)
          write(fvtPlun,'(1x,a,g12.5,a,g12.5)') 
     1      'invariant mass is ', m*1000, ' +/- ', sigm*1000
        end if

      end do

      return
      end
      subroutine doWu
C --------------------------------------------------------------------
      implicit none
      integer KRUN, KEVT
      integer NEVT
      include 'fvt.inc'

      integer i, i0, i1, new, status
      double precision x(3),Cx(3,3),d,sigd,m,sigm,ml(5)
      double precision chi2, p(3), Dp(3,3), Ep(3,3), pr
      integer fvtRead, fvtDo, fvtInvM
      external fvtRead, fvtDo, fvtInvM
      real PROB
      external PROB

      parameter (NEVT=16)
      integer runList(2,NEVT)
      integer trList(4,0:NEVT)
      integer pList(4,0:NEVT)
      double precision mass(4)
      data runList / 
     1  5166,000981,
     1  5825,000051,
     1  6972,004547,
     1  7422,002873,
     1  7496,003332,
     1  7681,002160,
     1  8126,000771,
     1  8163,002748,
     1  8334,007367,
     1  8335,001133,
     1  8537,004398,
     1  8623,000158,
     1  8833,004673,
     1  8855,003264,
     1  8865,006111,
     1  8960,000749
     1/
C -- tracks:
C -- event: Z0 -> l0 l2
C --              l0 -> l1 V+ V-
C -- V+, V-, l1, l2
      data trList /
     1  1,2,3,4,    
     1  4,3,2,1,
     1  3,4,1,2,
     1  1,4,2,3,
     1  3,2,4,1,
     1  1,2,4,3,
     1  4,1,2,3,
     1  3,4,1,2,
     1  3,2,1,4,
     1  2,4,1,3,
     1  2,3,4,1,
     1  3,4,2,1,
     1  4,2,3,1,
     1  3,2,4,1,
     1  4,1,3,2,
     1  1,2,4,3,
     1  4,1,2,3
     1/
C -- 1:electron, 2:muon, 3:tau, 4:pi
      data pList /
     1  4,4,4,4,    
     1  2,2,4,2,
     1  1,1,2,2,
     1  4,4,4,4,
     1  4,4,1,1,
     1  4,4,2,2,
     1  1,1,1,4,
     1  4,4,1,4,
     1  4,4,1,1,
     1  4,4,1,1,
     1  1,1,4,1,
     1  4,4,2,4,
     1  2,2,4,4,
     1  4,4,4,4,
     1  1,1,2,2,
     1  2,2,2,1,
     1  1,1,1,1
     1/
      data mass/
     1  0.00051099906d0,
     1  0.105658387d0,
     1  1.7841d0,
     1  0.1395675d0
     1/

      do i = 1, NEVT
        KRUN = runList(1,i)
        KEVT = runList(2,i)


C -- fit for common vertex of l+l-V=V-
        write(fvtPlun,'(1x,78(1H=))') 
        write(fvtPlun,'(/1x,a,i5,a,i5)') 
     1    'Run ', KRUN, ' Evt ', KEVT
        status = fvtRead(KRUN,KEVT)
        do i0 = 1,3
        do i1 = 1,3
          tCx(i0,i1,1) = tCx(i0,i1,1) * 10 000d0
        end do
        end do

        status = fvtDo(chi2, 4,trList(1,0))

C -- fit for vertex of V+V-
        status = fvtRead(KRUN,KEVT)
        do i0 = 1,3
        do i1 = 1,3
          tCx(i0,i1,1) = tCx(i0,i1,1) * 10 000d0
        end do
        end do
        write(fvtPlun,'(/1x,a,i5,a,i5)') 
     1    'common vertex for V+V- for Run ', KRUN, ' Evt ', KEVT

C -- fill tt, tCt from th, tCh for particles trList, pList
        call fvtHelix2P4(tt,tCt,
     1    2,trList(1,i),pList(1,i),th,tCh)
        status = fvtInvM(m,sigm,2,trList(1,i),tt,tCt)
        write(fvtPlun,'(1x,a,g12.5,a,g12.5)') 
     1    'invariant mass is ', m*1000, ' +/- ', sigm*1000
        status = fvtDo(chi2, 2,trList(1,i))
C -- save the V+V- vertex in x, Cx
        call fvCopy(x,tx,3,1)
        call fvCopy(Cx,tCx,3,3)

C -- calculate invariant mass of V+V-
        ml(trList(1,i)) = mass(pList(1,i))
        ml(trList(2,i)) = mass(pList(2,i))
        call fvInvM(m,sigM,
     1       2,trList(1,i),tw2pt,ml,x,Cx,tq,tCq,tGh)
        write(fvtPlun,'(1x,a,g12.5,a,g12.5)') 
     1    'invariant mass is ', m*1000, ' +/- ', sigm*1000

C -- sum up momenta of V+V- and put resultant track parameters into th/tGh
        new = tnt+1
        status = fvSumQ(p,Dp,Ep,
     1              new,tnt,trList,tw2pt,x,Cx,tq,tCq,th,tGh)

C -- fit vertex to l+l- pair
        write(fvtPlun,'(/1x,a,i5,a,i5)') 
     1    'common vertex for l(V)l for Run ', KRUN, ' Evt ', KEVT
        call fvtHelix2P4(tt,tCt,
     1    2,trList(3,i),pList(3,i),th,tCh)
        status = fvtInvM(m,sigm,2,trList(3,i),tt,tCt)
        write(fvtPlun,'(1x,a,g12.5,a,g12.5)') 
     1    'invariant mass is ', m*1000, ' +/- ', sigm*1000
        do i0 = 1,3
        do i1 = 1,3
          tCx(i0,i1,1) = tCx(i0,i1,1) * 10 000d0
        end do
        end do
        status = fvtDo(chi2, 2,trList(3,i))

        ml(trList(3,i)) = mass(pList(3,i))
        ml(trList(4,i)) = mass(pList(4,i))
        call fvInvM(m,sigM,
     1       2,trList(3,i),tw2pt,ml,tx(1,1),tCx(1,1,1),tq,tCq,tGh)
        write(fvtPlun,'(1x,a,g12.5,a,g12.5)') 
     1    'invariant mass is ', m*1000, ' +/- ', sigm*1000

        call fvtQM2P4(tt,tCt, 2,trList(3,i),pList(3,i),tq,tCq)
        status = fvtInvM(m,sigm,2,trList(3,i),tt,tCt)
        write(fvtPlun,'(1x,a,g12.5,a,g12.5)') 
     1    'invariant mass is ', m*1000, ' +/- ', sigm*1000

C -- calculate distance between vertex of V+V- and l+l-
        call fvDist(d, sigd, 
     1              x,Cx, tx(1,1), tCx(1,1,1))
        write(fvtPlun,'(/1x,a,g10.3,a,g10.3)') 
     1    'distance between V/ll vertices is ', d, ' +/- ', sigd

C -- fit distance between the two
        status = fvFitD(d,sigd,chi2,
     1                  d,tx(1,1),tCx(1,1,1),x,Cx,p,Dp,Ep)
        if (iand(status,1).NE.1) then
          print *, 'fit to tau did not work...'
        end if
        pr = PROB(REAL(chi2), 1)
        write(fvtPlun,'(1x,a,g10.3,a,g10.3,a,g10.3,a,g10.3)') 
     1    'fitted distance between V/ll vertices is ',
     1    d, ' +/- ', sigd, ' chi^2 ', chi2, ' prob. ', pr


      end do


      return
      end
      subroutine doOrs
C --------------------------------------------------------------------
      implicit none
      integer KRUN, KEVT
      integer NEVT
      include 'fvt.inc'
      integer i, i0, i1, status
      integer fvtRead, fvtDo, fvtInvM
      external fvtRead, fvtDo, fvtInvM
      double precision x(3),Cx(3,3),d,sigd,m,sigm
      double precision chi2
      parameter (NEVT=11)
      integer runList(2,NEVT)
      integer trList(4,0:NEVT)
      integer pList(4,0:NEVT)
      data runList / 
     1  4515,0750,
     1  5060,2137,
     1  5821,3675,
     1  5880,3296,
     1  7339,1745,
     1  7412,4101,
     1  7743,8370,
     1  8335,1133,
     1  8383,0966,
     1  8619,2198,
     1  8898,5940
     1/
C -- tracks:
C -- event: Z0 -> l0 l2
C --              l0 -> l1 V+ V-
C -- V+, V-, l1, l2
      data trList /
     1  1,2,3,4,    
     1  3,4,1,2,
     1  4,1,3,2,
     1  2,4,3,1,
     1  1,4,2,3,
     1  4,3,2,1,
     1  4,3,2,1,
     1  5,6,1,3,
     1  4,3,2,1,
     1  5,6,1,2,
     1  3,4,2,1,
     1  3,4,2,5
     1/
C -- 1:electron, 2:muon, 3:tau, 4:pi
      data pList /
     1  4,4,4,4,
     1  4,4,2,2,
     1  4,4,1,1,
     1  2,2,2,2,
     1  1,1,1,1,
     1  1,1,2,2,
     1  1,1,2,2,
     1  1,1,4,4,
     1  1,1,4,4,
     1  1,1,4,4,
     1  1,1,2,2,
     1  1,1,1,1
     1/

      do i = 1, NEVT
        KRUN = runList(1,i)
        KEVT = runList(2,i)


        write(fvtPlun,'(1x,78(1H=))') 
        write(fvtPlun,'(/1x,a,i5,a,i5)') 
     1    'Run ', KRUN, ' Evt ', KEVT
        status = fvtRead(KRUN,KEVT)
        do i0 = 1,3
        do i1 = 1,3
          tCx(i0,i1,1) = tCx(i0,i1,1) * 10 000d0
        end do
        end do

        status = fvtDo(chi2, 4,trList(1,0))

        status = fvtRead(KRUN,KEVT)
        do i0 = 1,3
        do i1 = 1,3
          tCx(i0,i1,1) = tCx(i0,i1,1) * 10 000d0
        end do
        end do
        write(fvtPlun,'(/1x,a,i5,a,i5)') 
     1    'common vertex for V+V- for Run ', KRUN, ' Evt ', KEVT

C -- fill tt, tCt from th, tCh for particles trList, pList
        call fvtHelix2P4(tt,tCt,
     1    2,trList(1,i),pList(1,i),th,tCh)
        status = fvtInvM(m,sigm,2,trList(1,i),tt,tCt)
        write(fvtPlun,'(1x,a,g12.5,a,g12.5)') 
     1    'invariant mass is ', m*1000, ' +/- ', sigm*1000
        status = fvtDo(chi2, 2,trList(1,i))
        call fvCopy(x,tx,3,1)
        call fvCopy(Cx,tCx,3,3)

        call fvtQM2P4(tt,tCt, 2,trList(1,i),pList(1,i),tq,tCq)
        status = fvtInvM(m,sigm, 2,trList(1,i),tt,tCt)
        write(fvtPlun,'(1x,a,g12.5,a,g12.5)') 
     1    'invariant mass is ', m*1000, ' +/- ', sigm*1000

        write(fvtPlun,'(/1x,a,i5,a,i5)') 
     1    'common vertex for l(V)l for Run ', KRUN, ' Evt ', KEVT
        call fvtHelix2P4(tt,tCt,
     1    2,trList(3,i),pList(1,i),th,tCh)
        status = fvtInvM(m,sigm,2,trList(3,i),tt,tCt)
        write(fvtPlun,'(1x,a,g12.5,a,g12.5)') 
     1    'invariant mass is ', m*1000, ' +/- ', sigm*1000
        status = fvtDo(chi2, 2,trList(3,i))
        call fvtQM2P4(tt,tCt, 2,trList(3,i),pList(1,i),tq,tCq)
        status = fvtInvM(m,sigm,2,trList(3,i),tt,tCt)
        write(fvtPlun,'(1x,a,g12.5,a,g12.5)') 
     1    'invariant mass is ', m*1000, ' +/- ', sigm*1000

C -- calculate distance
        call fvDist(d, sigd, 
     1              x, Cx, tx(1,1), tCx(1,1,1))
        write(fvtPlun,'(/1x,a,g10.3,a,g10.3)') 
     1    'distance between V/ll vertices is ', d, ' +/- ', sigd

      end do


      return
      end
      integer function fvtDo(chi2, nt,tList)
C-----------------------------------------------------------------------    
C -- do vertex fit.
C -- use tx(i,1) and tCx(i,j,1) as initial values
C -- and track parameters from th, tCh
C -- put results into tx, tCx and tq, tCq
C --
      implicit none
      double precision chi2
      integer nt, tList(nt)
      include 'fvt.inc'

      double precision x0(3), Cx0(3,3)

C -- ???? is not necessary
      call fvCopy(x0, tx,3,1)
      call fvCopy(Cx0, tCx,3,3)

      fvtDo = fvFit(tx,tCx,tq,tCq,tChi2,chi2,
     1    nt,tList,x0,Cx0,th,tGh)

      return
      end
      integer function fvtRead(KRUN, KEVT)
C-----------------------------------------------------------------------        
C -- read in track data, return
C --   tx, Cx: vertex position in carthesian coord. and cov. matrix
C --   th, Ch: track parameter in ALEPH comvention and cov. matrix
C --   tnt: number of tracks
C --   tGh: inverse of covariance matrix of track parameter
C -- in common /FVCOM/
C --
C -- returns ERROR if problem during inverting tCh occurred
C
      implicit none
      integer KRUN, KEVT
      include 'fvt.inc'
      character*80 fnam
      integer i0,i1,j, status

      fvtRead = NORMAL

      write(fnam, '(a,a,i5.5,a,i6.6,a)') 
     1  dataDir,'tr',KRUN,'e',KEVT,'.dat'
      open(unit=42, file = fnam,
     1 status = 'OLD')
      read(42,*) (tx(i0,1), i0=1,3)
      read(42,*) ((tCx(i0,i1,1), i0=1,3),i1=1,3)
      read(42,*) tw2pt
      read(42,*) tnt

      do j = 1, tnt
        read(42,*) (th(i0,j),i0=1,5)
        read(42,*) ((tCh(i0,i1,j),i0=1,5),i1=1,5)
      end do
      close(42)

C -- calculate inverse of covariance matrices for th
      do j = 1, tnt
CCCCCCCCCC        call fvCopy(tGh(1,1,j),tCh(1,1,j),5,5)
        status = fvLUinv(tGh(1,1,j), tCh(1,1,j),5)
CCCC        status = fvCalcG(tGh(1,1,j),tCh(1,1,j),5)
CCCC        if (status .EQ. INFO) then
CCCC          write(fvtPlun,'(1x,a)')  
CCCC     1      'editing of singular values during inv of Ch0'
CCCC          call fvAB(one,tGh(1,1,j),tCh(1,1,j),5,5,5)
CCCC          write(fvtPlun,'(/1x,5g10.3)') ((one(i0,i1), i0=1,5),i1=1,5)
CCCC        end if
        if (iand(status,1).NE.1) then
          fvtRead = status
        end if
      end do

      return
      end
      integer function fvtInvM(mi,sigmi, nt,tList,t,Ct)
C-----------------------------------------------------------------------        
C -- calculate invariant mass and error
C -- from 4-momenta t, Ct
C
      implicit none
      double precision mi, sigmi
      integer nt, tList(nt)
      double precision t(4,*), Ct(4,4,*)

      include 'fvt.inc'
      integer i, it, i0, i1
      double precision pi(4), Cpi(4,4)

      call fvZeroA(pi,4,1)
      call fvZeroA(Cpi,4,4)
      do i = 1, nt
        it = tList(i)
        do i0 = 1, 4
          pi(i0) = pi(i0) + t(i0,it)
          do i1 = 1, 4
            Cpi(i0,i1) = Cpi(i0,i1) + Ct(i0,i1,it)
          end do
        end do
      end do

      mi = sqrt(pi(4)**2-pi(1)**2-pi(2)**2-pi(3)**2)
      sigmi = 
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
      sigmi = sqrt(sigmi)/mi

      return
      end
      subroutine fvtHelix2P4(t,Ct, nt,tList,pList,p,Ch)
C-----------------------------------------------------------------------        
C -- calculate 4-momentum from track parameters
C -- using track/particle list tList, pList
C
      implicit none
      double precision t(4,*), Ct(4,4,*)
      integer nt,tList(nt),pList(nt)
      double precision p(5,*),Ch(5,5,*)

      include 'fvt.inc'

      integer i, it, ip
      double precision m
      double precision mass(4)
      data mass/
     1  0.00051099906d0,
     1  0.105658387d0,
     1  1.7841d0,
     1  0.1395675d0
     1/

      do i = 1, nt
        it = tList(i)
        call fvHelix2P4(t(1,it),Ct(1,1,it), 
     1             tw2pt,p(1,it),Ch(1,1,it))
        ip = pList(i)
        m = mass(ip)
        call fvEnergy(t(1,it),Ct(1,1,it), m,t(1,it),Ct(1,1,it))
      end do

      return
      end
      subroutine fvtQM2P4(tl,Ctl, nt,tList,pList,ql,Cql)
C-----------------------------------------------------------------------        
C -- calculate list of 4-momentum vectors 
C -- from mass and q-vectors {w, tl, psi} 
C -- using track/particle list tList, pList
C
      implicit none
      double precision tl(4,*), Ctl(4,4,*)
      integer nt,tList(nt),pList(nt)
      double precision ql(3,*),Cql(3,3,*)

      include 'fvt.inc'

      integer i, it, ip
      double precision m
      double precision mass(4)
      data mass/
     1  0.00051099906d0,
     1  0.105658387d0,
     1  1.7841d0,
     1  0.1395675d0
     1/

      do i = 1, nt
        it = tList(i)
        call fvQ2P4(tl(1,it),Ctl(1,1,it),ql(1,it),Cql(1,1,it),tw2pt)
        ip = pList(i)
        m = mass(ip)
        call fvEnergy(tl(1,it), Ctl(1,1,it), m)
      end do

      return
      end
      integer function fvtVert(v,Cv,q,Cq, v0,Cv0,q0,Cq0)
C-----------------------------------------------------------------------        
C --
      implicit none
      double precision v(3), Cv(3,3), q(3), Cq(3,3),
     1                 v0(3), Cv0(3,3), q0(3), Cq0(3,3)

      include 'fvt.inc'

      fvtVert = NORMAL

      return
      end
      integer function fvtTAV(x,Cx,p,Cp,chi2, x0,Cx0,h,Ch)
C-----------------------------------------------------------------------        
C -- track and vertex:
C -- do the refit of track with 
C -- helix parameters h and covariance matrix Ch,
C -- using point x0 and its covariance matrix Cx0
C -- h = {w, tl, psi0, d0, z0} and x0 = {x, y, z}
C -- return results in x, Cx, and p, Cp
C --
      implicit none
      double precision x(3), Cx(3,3), p(5), Cp(5,5), chi2
      double precision x0(3), Cx0(3,3), h(5), Ch(5,5)
      include 'fvt.inc'

      double precision v0(3), C0(3,3), q0(3), D0(3,3),
     1                 A(5,3), B(5,3), h0(5)
      double precision v(3), C(3,3), Gv(3,3), q(3), D(3,3), E(3,3)

      integer status, i

      fvtTAV = NORMAL

C -- calculate vertex and 3-momentum
C -- at approximate vertex position x0
C -- and the appropriate covariance matrices
      status = fvv0q0(v0, C0, q0, D0, x0, Cx0, h, Ch)
      if (iand(status,1).NE.1) then
        fvtTAV = status
        return
      end if

C -- calculate coefficients for measurement equation, A, B, h0
      status = fvABh0(A,B,h0, v0,q0)
      if (iand(status,1).NE.1) then
        fvtTAV = status
        return
      end if

C -- calculate filtered state vectors and covariance matrices
      
      status = fvFilter(v, C, Gv, q, D, E, chi2,
     1                  v0, C0, h, Ch, A, B, h0)
      if (iand(status,1).NE.1) then
        fvtTAV = status
        return
      end if

C -- re-calculate coefficients for measurement equation, A, B, h0
C -- and check on difference, ev. re-iterate (not done at the moment)
      status = fvABh0(A,B,h0, v,q)
      if (iand(status,1).NE.1) then
        fvtTAV = status
        return
      end if

C -- calculate estimated vertex and track parameters
      status = fvxp(x,Cx,p,Cp, v,q,C,D,E,A,B)

C -- some printout
      write(fvtPlun,'(1x,g10.3)')  chi2

      write(fvtPlun,'(1x,3(g10.3,a,g10.3))') (v0(i), ' +/-',
     1   sqrt(C0(i,i)) ,i=1,3)
      write(fvtPlun,'(1x,3(g10.3,a,g10.3))') (v(i), ' +/-',
     1   sqrt(C(i,i)) ,i=1,3)

CCC      write(fvtPlun,'(1x,3(g10.3,a,g10.3))') (x0(i), ' +/-',
CCC     1   sqrt(Cx0(i,i)) ,i=1,3)
CCC      write(fvtPlun,'(1x,3(g10.3,a,g10.3))') (x(i), ' +/-',
CCC     1   sqrt(Cx(i,i)) ,i=1,3)

      write(fvtPlun,'(1x,3(g10.3,a,g10.3))') (q0(i), ' +/-',
     1   sqrt(Ch(i,i)) ,i=1,3)
      write(fvtPlun,'(1x,3(g10.3,a,g10.3))') (q(i), ' +/-',
     1   sqrt(D(i,i)) ,i=1,3)

      write(fvtPlun,'(1x)')


      return
      end
      subroutine errInit()
      return
      end
      subroutine errSumm()
      return
      end
