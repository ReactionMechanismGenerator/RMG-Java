################################################################################
#
#	Makefile for DASSL
#
################################################################################

# The directory in which the source files can be found
$SOURCEDIR=.

# The directory in which to place temporary compiled files
$BUILDDIR=../../build/dassl

# The directory in which to place compiled executables and JAR files
$BINDIR=../../bin

# The Fortran 90 compiler to use and flags to use when compiling Fortran 90 
# code
F90=g95
F90FLAGS=

LDFLAGS=

OBJ=ddassl.o daux.o dlinpk.o call_dasslAUTO.o res.o getflux.o res_daepack.o reaction_flux.o

OBJ=$(patsubst %.o,$(BUILDDIR)/%.o,$(OBJ))

################################################################################

$(BINDIR)/dasslAUTO.exe: dirs $(OBJ)  
	$(F90) $(LDFLAGS) $(OBJ) -o $(BINDIR)/dasslAUTO.exe

dirs:
	mkdir -p $(BUILDDIR)
	mkdir -p $(BINDIR)

$(BUILDDIR)/%.o: %.f
	$(F90) $(F90FLAGS) -c $< -o $@

$(BUILDDIR)/%.o: %.f90
	$(F90) $(F90FLAGS) -c $< -o $@

clean:
	rm -rf $(BUILDDIR)
	rm -rf $(BINDIR)
