################################################################################
#
#	Makefile for DASPK
#
################################################################################

# The directory in which the source files can be found
$SOURCEDIR=.

# The directory in which to place temporary compiled files
$BUILDDIR=../../build/daspk

# The directory in which to place compiled executables and JAR files
$BINDIR=../../bin

# The Fortran 90 compiler to use and flags to use when compiling Fortran 90 
# code
F90=g95
F90FLAGS=

LDFLAGS=-lblas -llapack -ldaepack

OBJ=ddaspk.o dbanpre.o daux.o dilupre.o dlinpk.o drbdpre.o drbgpre.o dsensd.o dsparsk.o adf_dummy.o mpi_dummy.o call_daspkAUTO.o res.o getflux.o resad.o jac.o res_rpar.o calcdgdk.o reaction_flux.o

OBJ=$(patsubst %.o,$(BUILDDIR)/%.o,$(OBJ))

################################################################################

$(BINDIR)/daspkAUTO.exe: dirs $(OBJ)  
	$(F90) $(LDFLAGS) $(OBJ) -o $(BINDIR)/daspkAUTO.exe

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
