F90	= ifort
FLAGS	= -c -g
LFLAGS  = -ldaepack -llapack -lblas-3
FOBJECTS= ddaspk.o dbanpre.o daux.o dilupre.o dlinpk.o drbdpre.o drbgpre.o dsensd.o dsparsk.o adf_dummy.o mpi_dummy.o call_daspkAUTO.o res.o getflux.o resad.o jac.o res_rpar.o calcdgdk.o reaction_flux.o


dassl.exe:$(FOBJECTS)  
	$(F90) $(LFLAGS) $(FOBJECTS) -g -o daspkAUTO.exe
	cp daspkAUTO.exe $(RMG)/software/ODESolver
%.o:daspk/%.f
	$(F90) $(FLAGS) $<
%.o:%.f
	$(F90) $(FLAGS) $<
%.o:daspk/%.f90
	$(F90) $(FLAGS) $<
%.o:%.f90
	$(F90) $(FLAGS) $<
