# makefile for adiabatic reactor program

FC       = ifort
FFLAGS   = 
LINK     = ifort
LFLAGS   = -static-libcxa

OBJS     = driv.o reactor.o res.o parse.o readxml.o xmlparse.o
DASSLOBJ = ddassl.o daux.o dlinpk.o
CKOBJS   = cputim.o cklib.o

all: reactor.exe chem.exe

reactor.exe : $(OBJS) $(DASSLOBJ) $(CKOBJS)
	$(LINK) $(LFLAGS) -o reactor.exe $(OBJS) $(DASSLOBJ) $(CKOBJS)
	cp reactor.exe ../run

chem.exe : ckinterp.o
	$(LINK) $(LFLAGS) -o chem.exe ckinterp.o
	cp chem.exe ../run

readxml.o : readxml.f90
	$(FC) -c readxml.f90

xmlparse.o : xmlparse.f90
	$(FC) -c xmlparse.f90

.f.o:
	$(F77) -c $(FFLAGS) $<

clean :
	rm -f *.exe *.o

