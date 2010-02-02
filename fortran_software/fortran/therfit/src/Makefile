
F77    = ifort -g
FFLAGS = 
LINK   = ifort
LFLAGS = -static-libcxa

THOBJS = therfit.o intther.o covsrt.o gaussj.o mrqcof.o mrqmin.o

therfit.exe : $(THOBJS)
	$(LINK) $(LFLAGS) -o therfit.exe $(THOBJS)
	cp therfit.exe ../run

.f.o:
	$(F77) -c $(FFLAGS) $<

clean : 
	rm -f *.exe *.a *.o
