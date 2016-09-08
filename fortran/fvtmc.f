      subroutine fvtMCunify(p0, p)
C-----------------------------------------------------------------------        
C -- shift track parameters from p to p0 to a common vertex
C -- prelim.: by setting d0/z0 to zero
C
      implicit none
      double precision p0(5), p(5)
      integer i, it

      p0(1) = p(1)        
      p0(2) = p(2)        
      p0(3) = p(3)
      p0(4) = 1d-3        
      p0(5) = 1d-3        

      return
      end
      subroutine fvtMCsmear(p, p0,CpLam,CpU)
C-----------------------------------------------------------------------        
C -- smear track parameters p0 according to their error matrix
C -- denoted bu CpLam and CpU
C
      implicit none

      double precision p(5)
      double precision p0(5), CpLam(5), CpU(5,5)

      double precision s, pp(5)
      integer i, it, i0, i1
      real RG32
      external RG32

C -- get vector in diagonalized system, pp = U^T.p0
      do i0 = 1, 5
        s = 0d0
        do i1 = 1, 5
          s = s + CpU(i1,i0)*p0(i1)
        end do
        pp(i0) = s
      end do
C -- smear pp with diag{CpLam}
      do i0 = 1, 5
        pp(i0) = pp(i0)+CpLam(i0)*DBLE(RG32(i0))
      end do
C -- transform back p = U.pp
      do i0 = 1, 5
        s = 0d0
        do i1 = 1, 5
          s = s + CpU(i0,i1)*pp(i1)
        end do
        p(i0) = s
      end do

      return
      end
