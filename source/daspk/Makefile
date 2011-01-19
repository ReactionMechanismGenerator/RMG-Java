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

# The Fortran 90 compiler to use and flags to use when compiling Fortran code
# Call with 'make F90=g95' if you want to use g95 
# or 'make F90=gfortran' for GNU Fortran (recommended)
# Here we set the default
F90=g95
#F90=gfortran

ifeq ($(F90),g95)
F90FLAGS = -fmod=$(BUILDDIR) -ftrace=full   # used for dassl and daspk
endif ###### END OF g95 SETTINGS

ifeq ($(F90),gfortran)
# --fpe-trap=invalid,zero,underflow,overflow 
# The above flag causes problems at least with fame. Maybe a bug in gfortran? (Maybe a bug in fame.)
F90FLAGS = -ftrapv -fbounds-check -frange-check \
           -ggdb -J""$(BUILDDIR)"" -O3  -Wall -Wno-unused 
# if gfortran>4.3 then add -fbacktrace (it's not supported in earlier versions)
ifeq ($(shell gfortran --version 2>/dev/zero|grep -iqs '^GNU Fortran.* [4-9]\.[3-9]\.[0-9]' && echo "ok"), ok)
F90FLAGS += -fbacktrace
endif
F90FLAGS_NDEBUG = $(F90FLAGS) # used for dassl and daspk

# call as `make MACOS=true F90=gfortran` if you want MacOS X 10.6+ 64-bit intel core2 features
ifdef MACOS
F90FLAGS +=  -arch x86_64 -march=core2
F90_EXTRA_LDFLAGS +=  -framework vecLIB
endif

# Special flags for daspk files
# Don't check for bounds on arrays beacuse dlinpk breaks them on purpose.
# Don't warn about tab characters or uninitialized variables.
DASPK_FFLAGS = -fno-bounds-check -Wtabs -Wno-uninitialized
endif ###### END OF gfortran SETTINGS


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
# NB. they are also compiled with extra special flags
$(BUILDDIR)/%.o: private/%.f
	mkdir -p $(BUILDDIR)
	$(F90) $(F90FLAGS) $(DASPK_FFLAGS) -c $< -o $@

clean:
	rm -rf $(BUILDDIR)
	rm -rf $(BINDIR)/daspkAUTO.exe
