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
#F90=g95
#F90FLAGS=-fmod=$(BUILDDIR) -ftrace=full

F90=gfortran
F90FLAGS = -ffpe-trap=invalid,zero,overflow -ftrapv -fbounds-check -frange-check \
           -ggdb -J""$(BUILDDIR)"" -O3  -Wall -Wno-unused $(backtrace)
# if gfortran>4.3 then add -fbacktrace (it's not supported in earlier versions)
backtrace = $(shell gfortran --version 2>/dev/zero|grep -iqs '^GNU Fortran.* [4-9]\.[3-9]\.[0-9]' && echo "-fbacktrace")
F90FLAGS += -Wtabs # Don't warn about tab characters, although this breaks the Fortran standard.

LDFLAGS= $(F90_EXTRA_LDFLAGS) -lblas -llapack -ldaepack

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

# the files in the private folder are not part of the base RMG distribution
# because we are not licensed to distribute them
$(BUILDDIR)/%.o: private/%.f
	mkdir -p $(BUILDDIR)
	$(F90) $(F90FLAGS) -c $< -o $@
	
clean:
	rm -rf $(BUILDDIR)
	rm -rf $(BINDIR)/daspkAUTO.exe
