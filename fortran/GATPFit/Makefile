F77	= f77
F90	= g95
FLAGS	= -c -g -ftrace=full
FOBJECTS = daxpy.o   dgemm.o  dgerq2.o dlamch.o dlarf.o  dnrm2.o  dormr2.o dtrsm.o  lsfp_chem1.o xerbla.o dcopy.o  dgemv.o  dgerqf.o dlange.o dlarft.o dorm2r.o dormrq.o dtrsv.o  lsfp_nasa.o dgelq2.o dgeqr2.o dgglse.o dlapy2.o dlascl.o dorml2.o dscal.o  ieeeck.o lsfp_nist1.o dgelqf.o dgeqrf.o dggrqf.o dlarfb.o dlaset.o dormlq.o dtrmm.o  ilaenv.o lsfp.o dgels.o  dger.o  dlabad.o dlarfg.o dlassq.o dormqr.o dtrmv.o  lsame.o  lsfp_wilh1.o


GATPFit.exe:$(FOBJECTS)
	$(F90) $(FOBJECTS) -o GATPFit.exe
	cp GATPFit.exe $(RMG)/software/GATPFit
	
%.o:%.f
	$(F90) $(FLAGS) $<

clean:
	rm -f *.o *.exe
