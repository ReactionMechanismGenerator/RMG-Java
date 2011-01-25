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

# The Fortran 90 compiler to use and flags to use when compiling Fortran code
# Call with 'make F90=g95' if you want to use g95 
# or 'make F90=gfortran' for GNU Fortran (recommended)
# Here we set the default
#F90=g95
F90=gfortran

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

# Special flags for dassl d*.f files
# Don't check for bounds on arrays beacuse dlinpk breaks them on purpose.
# Don't warn about tab characters
DASSL_FFLAGS = -fno-bounds-check -Wtabs
endif ###### END OF gfortran SETTINGS


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
