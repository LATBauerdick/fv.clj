C----------------------------------------------------------------------
C
C FVA, the ALPHA interface to FV, a package to fit vertices and find 
C secondary vertices
C
C author: Lothar A.T. Bauerdick, CERN/PPE
C date:   Jan. 1991
C
C----------------------------------------------------------------------
      subroutine fva
      end
      subroutine fvaIni()
C-----------------------------------------------------------------------
C -- initialize fv package
C -- for VAX processors: initialize error handler
C -- default: switch on multiple scattering between
C -- innermost coordinate of track fit and origin
C
      implicit none

      logical doMS
      common /FVACOM/ doMS

CVAX  call fvVaxEI()
      doMS = .true.
      
      end
      integer function fvaSMS(onOff)
C-----------------------------------------------------------------------
C -- switch on/off treatment of multiple scattering between
C -- last coordinate and beam spot
C
      implicit none
      integer onOff
      integer NORMAL, ERROR
      parameter (NORMAL=1, ERROR=2)

      logical doMS
      common /FVACOM/ doMS

      if (onOff .EQ. 0) then
         doMS = .false.
         call UFSWIM(0)
      else
         doMS = .true.
         call UFSWIM(1)
      end if
      fvaSMS = NORMAL
      return 
      end
      integer function fvaRFH(nt,tList)
C-----------------------------------------------------------------------   
C -- refit helices for nt tracks in tList (ALPHA track numbers)
C
      implicit none
      integer nt, tList(nt)

      integer NORMAL, ERROR
      parameter (NORMAL=1, ERROR=2)

      integer i, status
      integer fvaRefT
      external fvaRefT

      fvaRFH = NORMAL
      do i = 1, nt
        status = fvaRefT(tList(i))
        if (iand(status,1).NE.1) fvaRFH = status
      end do

      return
      end
      integer function fvaGHG(h,Gh, nt,tList)
C----------------------------------------------------------------------- 
C -- fill arrays h, Gh with track parameters and inverted
C -- covariance matrices for nt tracks in tList
C -- track list tList contains ALPHHA track numbers
C
      implicit none
      integer nt, tList(nt)
      double precision h(5,nt), Gh(5,5,nt)

      double precision hh(3),Ch(5,5)
      integer i, it, fvaGetH, fvLUinv, status, fvSVDinv
      external fvaGetH, fvLUinv, fvSVDinv

      integer NORMAL
      parameter (NORMAL=1)

      fvaGHG = NORMAL

      do i = 1, nt
        it = tList(i)
        status = fvaGetH(h(1,i),Ch,hh, it)
        if (iand(status,1).NE.1) then
C -- error in getting track parameters
          fvaGHG = status
        end if
CCCCCCC        status = fvLUinv(Gh(1,1,i), Ch,5)
        call fvCopy(Gh(1,1,i), Ch,5,5)
        status = fvSVDinv(Gh(1,1,i), 5)
        call fvTinv(Gh(1,1,i),Ch,5)
        if (iand(status,1).NE.1) then
          fvaGHG = status
        end if
        
      end do

      return
      end
      integer function fvaGHC(h,Ch, nt,tList)
C-----------------------------------------------------------------------       C -- fill arrays h, Ch with track parameters and covariance matrices
C -- for nt tracks in tList
C -- track list tList contains ALPHA track numbers
C
      implicit none
      double precision h(5,*), Ch(5,5,*)
      integer nt, tList(nt)

      double precision hh(3)
      integer i, it, fvaGetH, status
      external fvaGetH

      integer NORMAL
      parameter (NORMAL=1)

      fvaGHC = NORMAL

      do i = 1, nt
        it = tList(i)
        status = fvaGetH(h(1,i),Ch(1,1,i),hh, it)
        if (iand(status,1).NE.1) then
C -- error in getting track parameters

        end if
      end do

      return
      end
      integer function fvaGetH(p,Cp,pp, ita)
C-----------------------------------------------------------------------
C -- get track parameter p = {w, tl, psi0, d0, z0}, 
C -- its covariance matrix Cp = cov{p,p} and 
C -- additional parameter pp = {chi2, ndof, ier}
C -- for ALPHA track # ita
C -- 
      implicit none
      double precision p(5), Cp(5,5), pp(3)
      integer ita

      include '/aleph/include/alpha111/qdecl.inc'
      include '/aleph/include/alpha111/qcde.inc'

      integer NORMAL, ERROR
      parameter (NORMAL=1, ERROR=2)

      logical doSMS

      logical doMS
      common /FVACOM/ doMS

      integer kFRFT, ioff, nTr, iErr
      integer i0,i1,emOff, x, y, it
      integer nITC, nVDET
      real triag(15)

      integer NLINK
      external NLINK
C
C -- offset for track i in 'FRFT'
      integer KKFRFT
      KKFRFT(i0) = kFRFT+2+(i0-1)*(IW(kFRFT+1))

      include '/aleph/include/alpha111/qmacro.inc'

      fvaGetH = NORMAL

C -- if do multiple scattering try to find FRFT#2
C -- if not there do simple minded multiple scattering
C -- (later: re-fit tracks w/ swim to origin)
      doSMS = .false.
      if (doMS) then
         kFRFT = NLINK('FRFT',2)
         if (kFRFT .eq. 0) then
            kFRFT = NLINK('FRFT',0)
            doSMS = .true.
         end if
      else
         kFRFT = NLINK('FRFT',0)
      end if
C -- is bank there?
      if (kFRFT .eq. 0) then
        fvaGetH = ERROR
        return
      end if

C -- calculate FRFT row from ALPHA track number
      it = KTN(ita)
C -- number of tracks?
      nTr = IW(kFRFT+2)
      if (it .gt. nTr) then
        fvaGetH = ERROR
        return
      end if

C -- get parameters for track it
      ioff = KKFRFT(it)

C -- get error, chi2 and # degree of freedom from track fit
      iErr = IW(ioff+30)
      iErr = iErr - 10*(iErr/10)
      pp(3) = iErr
      if (it .gt. nTr) then
        fvaGetH = ERROR
        return
      end if
      pp(1) = RW(ioff + 28)
      pp(2) = IW(ioff + 29)

C -- get p = {w, tl, psi0, d0, z0}
      do  i0 = 1, 5
        p(i0) = RW(ioff + i0)
      enddo

C -- do simple minded multiple scattering between innermost coordinate
C -- and the interaction point region if flag set
      if (doSMS) then
         nITC = KFRTNI(ita)
         nVDET = KFRTNV(ita)
         call UMSERR(1,nITC,nVDET,RW(ioff+1),RW(ioff+7),
     &            triag,iErr)
         if (iErr .ne. 0) then
           fvaGetH = ERROR
           return
         end if
      end if

C -- get Cp = cov{w, tl, psi0, d0, z0}
C -- covariance matrix stored in order
C --    1  2  4  7 11 16
C --       3  5  8 12 17
C --          6  9 13 18
C --            10 14 19
C --               15 20
C --                  21
      do i0 = 1, 5
      do i1 = 1, 5
        x = max(i0,i1)-1
        y = min(i0,i1)-1
        emOff = (x*(x+1))/2+y
        Cp(i0,i1) = triag(1+emOff)
      end do
      end do

      return
      end
      integer function fvaRefT(ita)
C-----------------------------------------------------------------------
C -- refit ALPHA track number ita 
C -- calls UFTTRK/UFTKAL from ALEPHLIB
C -- 
      implicit none
      integer ita

      integer NORMAL, ERROR
      parameter (NORMAL=1, ERROR=2)

      include '/aleph/include/alpha111/qdecl.inc'
      include '/aleph/include/alpha111/qcde.inc'

      logical first
      data first/.true./
      integer nFRTL, nFTCL, nFICL, nFVCL
      SAVE first, nFRTL, nFTCL, nFICL, nFVCL

      integer j, kFRFT, kFRTL, kFTCL, kFICL, kFVCL
      integer NTPC, NITC, NVDT
      integer NPTPC, NPITC, NPVDT

      real VV0(6),COV(21),CHI2
      integer NDEG

      integer iFRFT
      integer IERR, k

      integer NAMIND, NLINK
      external NAMIND, NLINK

      include '/aleph/include/alpha111/qmacro.inc'

      fvaRefT = NORMAL

      if (first) then
        first = .false.
        nFRTL = NAMIND('FRTL')
        nFTCL = NAMIND('FTCL')
        nFICL = NAMIND('FICL')
        nFVCL = NAMIND('FVCL')
      end if

      kFRFT = NLINK('FRFT',0)
      kFRTL = IW(nFRTL)
      kFTCL = IW(nFTCL)
      kFICL = IW(nFICL)
      kFVCL = IW(nFVCL)

      if (kFRFT.EQ.0 .OR.
     1    kFRTL.EQ.0 .OR.
     1    kFTCL.EQ.0 .OR.
     1    kFICL.EQ.0 .OR.
     1    kFVCL.EQ.0) then
        return
      end if

      J = KTN(ita)

      NTPC = ITABL(kFRTL,J,JFRTNT)
      NITC = ITABL(kFRTL,J,JFRTNI)
      NVDT = ITABL(kFRTL,J,JFRTNV)

      if (NVDT+NTPC .EQ. 0) then
        return
      end if
      if (NTPC .EQ. 0) then
      end if

      NPTPC = ITABL(kFRTL,J,JFRTIT)
      NPITC = ITABL(kFRTL,J,JFRTII)
      NPVDT = ITABL(kFRTL,J,JFRTIV)

      call UFTTRK(
     1  QMFLD, 
     1  RW(KROW(kFRFT,J)+JFRFIR),
     1  RTABL(kFRFT,J,JFRFC2),
     1  NTPC, NITC, NVDT,
     1  IW(kFTCL+LMHLEN+NPTPC+1),
     1  IW(kFICL+LMHLEN+NPITC+1),
     1  IW(kFVCL+LMHLEN+NPVDT+1),
     1  VV0, COV, CHI2, NDEG)

C -- fill FRFT bank
      iFRFT = KROW(kFRFT,J)
      if (CHI2 .GT. 1e10) then
        IERR = 4
        fvaRefT = ERROR
      else
        IERR = 0
        IW(iFRFT+JFRFDF) = NDEG
        IW(iFRFT+JFRFNO) = 90 + IERR
      end if
      if (IERR .LT. 4) then
        do k=1, 6
          RW(iFRFT+k) = VV0(k)
        end do
        do k=1, 21
          RW(iFRFT+6+k) = COV(k)
        end do
        RW(iFRFT+JFRFC2) = CHI2
      end if

      return
      end
