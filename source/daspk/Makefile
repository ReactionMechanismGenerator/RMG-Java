################################################################################
#
#	Makefile for DASPK
#
################################################################################

# The directory in which the source files can be found
SOURCEDIR=.

# The directory in which to place temporary compiled files
BUILDDIR=../../build/daspk

# The directory in which to place compiled executables and JAR files
BINDIR=../../bin

# The Fortran 90 compiler to use and flags to use when compiling Fortran 90 
# code
F90=g95
F90FLAGS=-fmod=$(BUILDDIR) -ftrace=full

LDFLAGS=-lblas -llapack -ldaepack

OBJ0=ddaspk.o dbanpre.o daux.o dilupre.o dlinpk.o drbdpre.o drbgpre.o dsensd.o dsparsk.o adf_dummy.o mpi_dummy.o call_daspkAUTO.o res.o getflux.o resad.o jac.o res_rpar.o calcdgdk.o reaction_flux.o

OBJ=$(patsubst %.o,$(BUILDDIR)/%.o,$(OBJ0))

################################################################################

$(BINDIR)/daspkAUTO.exe: $(OBJ)  
	mkdir -p $(BINDIR)
	$(F90) $(LDFLAGS) $(OBJ) -o $(BINDIR)/daspkAUTO.exe

$(BUILDDIR)/%.o: %.f
	mkdir -p $(BUILDDIR)
	$(F90) $(F90FLAGS) -c $< -o $@

$(BUILDDIR)/%.o: %.f90
	mkdir -p $(BUILDDIR)
	$(F90) $(F90FLAGS) -c $< -o $@

clean:
	rm -rf $(BUILDDIR)
	rm -rf $(BINDIR)/daspkAUTO.exe
