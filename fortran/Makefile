NAME    = FV
CFLAGS  = -w
FFILES  = fvt.f fvtmc.f prob.f erf64.f rg32.f rn32.f # doHad.f doMc.f fvtMCdiag.f # fva.f
FFLAGS  = -O
FC = gfortran
LDFLAGS =
ALIBS    = \
	  	/usr/local/lib/alpha.o \
	  	-lalpha \
		-lmini \
		-lalephlib \
		-lbos77

#

#LATB LIBS    = -L/cern/pro/lib \
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

fv.o : fv.f fv.inc

fva.o : fva.f

fvt.o : fvt.f fvt.inc

fvtmc.o : fvtmc.f fvt.inc

