      subroutine fvtMCdiag(CpLam,CpU, Cp)
C-----------------------------------------------------------------------        
C -- diagonalize covariance matrix Cp and return
C -- matrix CpU = {eigen vectors} and CpLam = sqrt(eigen values}
C --
      implicit none
      double precision CpLam(5),CpU(5,5), Cp(5,5)
      real AR(5,5), WR(5), ZR(5,5), WORK(5)
      integer i0, i1, ierr

      do i0 = 1, 5
      do i1 = 1, 5
        AR(i0,i1) = REAL(Cp(i0,i1))
      end do
      end do

      call EISRS1(5, 5, AR, WR, ZR, IERR, WORK)

      do i0 = 1, 5
      if (WR(i0) .GE. 0.) then
        CpLam(i0) = DBLE(sqrt(WR(i0)))
      else
        CpLam(i0) = 0d0
        call fvError('sqrt of negative number in fvtMCdiag')
      end if
      do i1 = 1, 5
        CpU(i0,i1) = DBLE(ZR(i0,i1))
      end do
      end do

      return
      end
