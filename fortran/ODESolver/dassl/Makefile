F77	= f77
FLAGS	= -c -g
LFLAGS  = 
FOBJECTS= ddassl.o daux.o dlinpk.o call_dassl.o res.o getflux.o res_daepack.o reaction_flux.o


dassl.exe:$(FOBJECTS)  
	$(F77) $(LFLAGS) $(FOBJECTS) -g -o dassl.exe
	cp dassl.exe $(RMG)/software/ODESolver
%.o:%.f
	$(F77) $(FLAGS) $<

