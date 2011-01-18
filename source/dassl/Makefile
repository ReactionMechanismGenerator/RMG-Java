################################################################################
#
#	Makefile for DASSL
#
################################################################################

# The directory in which the source files can be found
SOURCEDIR=.

# The directory in which to place temporary compiled files
BUILDDIR=../../build/dassl

# The directory in which to place compiled executables and JAR files
BINDIR=../../bin

# The Fortran 90 compiler to use and flags to use when compiling Fortran 90 
# code
#F90=g95
#F90FLAGS=-fmod=$(BUILDDIR) -ftrace=full

F90=gfortran
# --fpe-trap=underflow,overflow causes problems at least with fame. Maybe a bug in gcc? (Maybe a bug in fame.)
F90FLAGS = -ffpe-trap=invalid,zero -ftrapv -frange-check \
           -ggdb -J""$(BUILDDIR)"" -O3  -Wall -Wno-unused $(backtrace)
# if gfortran>4.3 then add -fbacktrace (it's not supported in earlier versions)
backtrace = $(shell gfortran --version 2>/dev/zero|grep -iqs '^GNU Fortran.* [4-9]\.[3-9]\.[0-9]' && echo "-fbacktrace")

# Special flags for dassl d*.f files
# Don't check for bounds on arrays beacuse dlinpk breaks them on purpose.
# Don't warn about tab characters
DASSL_FFLAGS = -fno-bounds-check -Wtabs

LDFLAGS= $(F90_EXTRA_LDFLAGS)

OBJ0=ddassl.o daux.o dlinpk.o call_dasslAUTO.o res.o getflux.o res_daepack.o reaction_flux.o

OBJ=$(patsubst %.o,$(BUILDDIR)/%.o,$(OBJ0))

################################################################################

$(BINDIR)/dasslAUTO.exe: $(OBJ)  
	mkdir -p $(BINDIR)
	$(F90) $(LDFLAGS) $(OBJ) -o $(BINDIR)/dasslAUTO.exe

# only .f files that begin with a d (eg. daux.f, dassl.f, dlinpk.f)
$(BUILDDIR)/d%.o: d%.f
	mkdir -p $(BUILDDIR)
	$(F90) $(F90FLAGS) $(DASSL_FFLAGS) -c $< -o $@

$(BUILDDIR)/%.o: %.f
	mkdir -p $(BUILDDIR)
	$(F90) $(F90FLAGS) -c $< -o $@

$(BUILDDIR)/%.o: %.f90
	mkdir -p $(BUILDDIR)
	$(F90) $(F90FLAGS) -c $< -o $@

clean:
	rm -rf $(BUILDDIR)
	rm -rf $(BINDIR)/dasslAUTO.exe
