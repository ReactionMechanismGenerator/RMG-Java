F77	= f77
FLAGS	= -c -g
LFLAGS  = -ldaepack -llapack -lblas-3
FOBJECTS= ddaspk.o dbanpre.o daux.o dilupre.o dlinpk.o drbdpre.o drbgpre.o dsensd.o dsparsk.o adf_dummy.o mpi_dummy.o call_daspk.o res.o getflux.o resad.o jac.o res_rpar.o calcdgdk.o reaction_flux.o


dassl.exe:$(FOBJECTS)  
	$(F77) $(LFLAGS) $(FOBJECTS) -g -o daspk.exe
	cp daspk.exe $(RMG)/software/ODESolver
%.o:daspk/%.f
	$(F77) $(FLAGS) $<
%.o:%.f
	$(F77) $(FLAGS) $<

