
F77    = ifort
FFLAGS = -O0 -lowercase
LINK   = ifort
LFLAGS = -O0 -lowercase -static-libcxa

CDOBJS = cdaux.o cddiss.o cdfouts.o cdkouts.o cdsets.o cdsubs.o \
	 cdbits.o cdgets.o cdstab.o chemdis2.o gammln.o\
	 cdchem.o cdkfits.o cdrecur.o cdstr.o genRhos.o straightBS.o \
	 PrintBarkerFiles.o PrintSteadyInput.o chebyshev.o

chemdis.exe : $(CDOBJS)
	$(LINK) $(LFLAGS) -o chemdis.exe $(CDOBJS)
	cp chemdis.exe ../run

.f.o:
	$(F77) -c $(FFLAGS) $<

clean : 
	rm -f *.exe *.a *.o
