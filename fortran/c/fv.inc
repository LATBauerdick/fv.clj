
      integer floatSize
      parameter (floatSize = 2)
      
      integer NORMAL, ERROR, INFO
      parameter (NORMAL=1, ERROR = 2, INFO = 3)
      double precision RR0
      parameter (RR0 = 0d0)

      integer dh,dv,dq
      parameter (dh=5,dv=3,dq=3)
      character*132 text
      logical print
      integer plun
      common /FVCOM/ 
     1                  print,
     1                  plun,
     1                  text


