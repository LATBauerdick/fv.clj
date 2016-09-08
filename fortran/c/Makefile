NAME    = FV
CFLAGS  = -w -I/sw/include/ -m32
FFILES  =  # doHad.f doMc.f fvtMCdiag.f # fva.f
CFILES  = fvt.c fvtmc.c prob.c erf64.c rg32.c rn32.c 
FFLAGS  = -O
FC = gfortran
LDFLAGS =  -m32
ALIBS    = \
	  	/usr/local/lib/alpha.o \
	  	-lalpha \
		-lmini \
		-lalephlib \
		-lbos77

#

LIBS    = -L/cern/pro/lib -L/sw/lib/ -lf2c -lm\
#		-lmathlib \
#		-lpacklib \
#		-lkernlib
OFILE_DIR = obj

# Rules...

SRCFILES = $(CFILES) $(FFILES)
OBJFILES = $(CFILES:.c=.o) $(FFILES:.f=.o) 

###$(FC) $(FFLAGS) 

$(NAME): $(OBJFILES) fv.o 
	$(FC) -o $@ fv.o $(OBJFILES) $(LIBS) $(LDFLAGS)

fve : fve.o fv.o fva.o 
	$(FC) -o $@ fve.o fv.o fva.o $(ALIBS) $(LIBS) $(LDFLAGS)

ks : ks.o fvks.o badTrack.o fv.o fva.o 
	$(FC) -o $@ ks.o fv.o $(ALIBS) $(LIBS) $(LDFLAGS)

fve.o : fve.f

ks.o : ks.f

fv.o : fv.c

fva.o : fva.f

fvt.o : fvt.c 

fvtmc.o : fvtmc.c

