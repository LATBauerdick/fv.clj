      integer NORMAL, ERROR, INFO
      parameter (NORMAL=1, ERROR = 2, INFO = 3)
      integer MAXTR
      parameter (MAXTR=40)

      character*(*) dataDir
      parameter (dataDir = 'dat/')
CCC      parameter (dataDir = 'disk$user:[bauerdick.fv.dat]')
      
      double precision tw2pt
      
      double precision th(5,MAXTR), tCh(5,5,MAXTR), 
     1                 tx(3,MAXTR), tCx(3,3,MAXTR),
     1                 tGh(5,5,MAXTR),
     1                 tt(4,MAXTR), tCt(4,4,MAXTR),
     1                 tq(3,MAXTR), tCq(3,3,MAXTR), 
     1                 tChi2(MAXTR)
      integer tnt
      logical fvtPrint
      integer fvtPlun
      common /FVTCOM/
     1              tw2pt,
     1              th, tCh, tx, tCx, 
     1              tGh,
     1              tt, tCt, 
     1              tq, tCq, 
     1              tChi2,
     1              tnt,
     1              fvtPrint,
     1              fvtPlun
      integer fvSumQ, fvFitD, fvFit
      external fvSumQ, fvFitD, fvFit
      double precision fvProb
      external fvProb
      integer fvCalcG, fvLUinv
      external fvCalcG, fvLUinv
      integer fvv0q0, fvxp
      external fvv0q0, fvxp
      integer fvABh0, fvFilter
      external fvABh0, fvFilter
