F90	= g95
FLAGS	= -c -g
LFLAGS  = 
FOBJECTS= ddassl.o daux.o dlinpk.o call_dasslAUTO.o res.o getflux.o res_daepack.o reaction_flux.o


dasslAUTO.exe:$(FOBJECTS)  
	$(F90) $(LFLAGS) $(FOBJECTS) -g -o dasslAUTO.exe
	cp dasslAUTO.exe $(RMG)/software/ODESolver
%.o:%.f
	$(F90) $(FLAGS) $<
%.o:%.f90
	$(F90) $(FLAGS) $<

clean:
	rm -f *.o *.exe
