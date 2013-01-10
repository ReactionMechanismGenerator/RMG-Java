################################################################################
#
#	Makefile for RMG
#
################################################################################

# The directory in which the source files can be found
SOURCEDIR=$(CURDIR)/source

# The directory in which to place temporary compiled files
BUILDDIR=$(CURDIR)/build

# The directory in which to place compiled executables and JAR files
BINDIR=$(CURDIR)/bin

# The directory in which to run RMG (used for the test)
RUNDIR=$(CURDIR)/run

# The Fortran 90 compiler to use and flags to use when compiling Fortran code
# Call with 'make F90=g95' if you want to use g95 
# or 'make F90=gfortran' for GNU Fortran (recommended)
# Here we set the default
F90=gfortran

# Some tests
uname_S := $(shell sh -c 'uname -s 2>/dev/null || echo not')
uname_M := $(shell sh -c 'uname -m 2>/dev/null || echo not')
f90version := $(shell $(F90) --version 2>/dev/zero)
gfortranversion := $(shell expr "$(f90version)" : 'GNU Fortran.*\([0-9]\.[0-9]\.[0-9]\)' )

# Detect 64 bit MacOS X for some optimizations
ifeq ($(uname_S),Darwin)
  ifeq ($(uname_M),x86_64)
    MACOS = True
  endif
endif

# If running G95
ifeq ($(shell expr "$(f90version)" : 'G95'), 3)
  F90FLAGS = -fbounds-check -ftrace=full -fmod=$(BUILDDIR) -Wall -O3
  F90FLAGS_NDEBUG = -fmod=$(BUILDDIR) -ftrace=full   # used for dassl and daspk
endif ###### END OF G95 SETTINGS

# If running gfortran
ifeq ($(shell expr "$(f90version)" : 'GNU Fortran'), 11)
  # --fpe-trap=invalid,zero,underflow,overflow 
  # The above flag causes problems at least with fame. Maybe a bug in gfortran? (Maybe a bug in fame.)
  F90FLAGS = -ftrapv -fbounds-check -frange-check \
             -ggdb -J""$(BUILDDIR)"" -O3  -Wall -Wno-unused

  # if gfortran >= 4.3 then add -fbacktrace (it's not supported in earlier versions)
  ifeq ($(shell expr "$(f90version)" : '[4-9]\.[3-9]'),3)
    F90FLAGS += -fbacktrace
  endif

  # Call as `make MACOS=true F90=gfortran` if you want MacOS X 10.6+ 64-bit intel core2 features
  ifdef MACOS
    F90FLAGS +=  -arch x86_64 -march=core2
    F90_EXTRA_LDFLAGS +=  -framework vecLIB
  endif

  F90FLAGS_NDEBUG = $(F90FLAGS) # used for dassl and daspk
endif ###### END OF gfortran SETTINGS

# these are added to the LDFLAGS of the subsidiary makefiles
#F90_EXTRA_LDFLAGS = -L/home/local/lib -lg2c  # required for Monch

################################################################################

# Default is to build those that come complete with the RMG distribution
base: dirs fame frankie GATPFit dassl RMG

# Make just the Fortran dependencies (i.e. not the Java)
fortran: dirs fame frankie GATPFit dassl

# Make all the Fortran (including DASPK)
all_fortran: dirs fame frankie GATPFit dassl daspk

# You can also build everything
all: dirs fame frankie GATPFit dassl daspk RMG inchi symmetry

RMG: dirs 
	mkdir -p $(BUILDDIR)/RMG
	ant clean # we don't trust ant to spot what needs doing!
	#ant compile
	ant jar

fame: dirs
	make -C $(SOURCEDIR)/fame SOURCEDIR=$(SOURCEDIR)/fame BUILDDIR=$(BUILDDIR)/fame BINDIR=$(BINDIR) F90=$(F90) F90FLAGS="$(F90FLAGS)" F90_EXTRA_LDFLAGS="$(F90_EXTRA_LDFLAGS)"

frankie: dirs
	make -C $(SOURCEDIR)/frankie SOURCEDIR=$(SOURCEDIR)/frankie BUILDDIR=$(BUILDDIR)/frankie BINDIR=$(BINDIR) F90=$(F90) F90FLAGS="$(F90FLAGS)" F90_EXTRA_LDFLAGS="$(F90_EXTRA_LDFLAGS)"

dassl: dirs
	make -C $(SOURCEDIR)/dassl SOURCEDIR=$(SOURCEDIR)/dassl BUILDDIR=$(BUILDDIR)/dassl BINDIR=$(BINDIR) F90=$(F90) F90FLAGS="$(F90FLAGS_NDEBUG)" F90_EXTRA_LDFLAGS="$(F90_EXTRA_LDFLAGS)"

daspk: dirs
	make -C $(SOURCEDIR)/daspk SOURCEDIR=$(SOURCEDIR)/daspk BUILDDIR=$(BUILDDIR)/daspk BINDIR=$(BINDIR) F90=$(F90) F90FLAGS="$(F90FLAGS_NDEBUG)" F90_EXTRA_LDFLAGS="$(F90_EXTRA_LDFLAGS)"

GATPFit: dirs
	make -C $(SOURCEDIR)/GATPFit SOURCEDIR=$(SOURCEDIR)/GATPFit BUILDDIR=$(BUILDDIR)/GATPFit BINDIR=$(BINDIR) F90=$(F90) F90FLAGS="$(F90FLAGS)" F90_EXTRA_LDFLAGS="$(F90_EXTRA_LDFLAGS)"

dirs:
	mkdir -p $(BUILDDIR)
	mkdir -p $(BINDIR)


# GET AND BUILD INCHI SUPPORT
inchi: $(BINDIR)/cInChI-1

$(BINDIR)/cInChI-1: $(BUILDDIR)/InChI-1-software-1-02-beta
	make -C $(BUILDDIR)/InChI-1-software-1-02-beta/INCHI_API/gcc_makefile
	cp $(BUILDDIR)/InChI-1-software-1-02-beta/INCHI_API/gcc_makefile/cInChI-1 $(BINDIR)/

$(BUILDDIR)/InChI-1-software-1-02-beta: $(BUILDDIR)/inchi102b.zip
	cd $(BUILDDIR); unzip inchi102b.zip

$(BUILDDIR)/inchi102b.zip:
	cd $(BUILDDIR); wget http://old.iupac.org/inchi/download/inchi102b.zip
# END OF INCHI SUPPORT

# GET AND BUILD SYMMETRY
symmetry: $(BINDIR)/SYMMETRY.EXE

$(BINDIR)/SYMMETRY.EXE: $(BUILDDIR)/symmetry/symmetry.c
	cd $(BUILDDIR)/symmetry/; cc -o $(BINDIR)/SYMMETRY.EXE -O3 -ansi -Wall symmetry.c -lm

$(BUILDDIR)/symmetry/symmetry.c: $(BUILDDIR)/symmetry.zip
	cd $(BUILDDIR); unzip symmetry.zip -d symmetry

$(BUILDDIR)/symmetry.zip:
	cd $(BUILDDIR); wget http://www.cobalt.chem.ucalgary.ca/ps/symmetry/symmetry.zip
# END OF SYMMETRY SUPPORT

clean:
	make -C $(SOURCEDIR)/fame clean SOURCEDIR=$(SOURCEDIR)/fame BUILDDIR=$(BUILDDIR)/fame BINDIR=$(BINDIR) 
	make -C $(SOURCEDIR)/frankie clean SOURCEDIR=$(SOURCEDIR)/frankie BUILDDIR=$(BUILDDIR)/frankie BINDIR=$(BINDIR) 
	make -C $(SOURCEDIR)/dassl clean SOURCEDIR=$(SOURCEDIR)/dassl BUILDDIR=$(BUILDDIR)/dassl BINDIR=$(BINDIR) 
	make -C $(SOURCEDIR)/daspk clean SOURCEDIR=$(SOURCEDIR)/daspk BUILDDIR=$(BUILDDIR)/daspk BINDIR=$(BINDIR) 
	make -C $(SOURCEDIR)/GATPFit clean SOURCEDIR=$(SOURCEDIR)/GATPFit BUILDDIR=$(BUILDDIR)/GATPFit BINDIR=$(BINDIR) 
	ant clean
	rm -rf $(BUILDDIR)
	rm -rf $(BINDIR)
	rm -rf $(RUNDIR)

# Run a test case
test:
	mkdir -p $(RUNDIR)
	cp examples/RMG/hexadiene/condition.txt $(RUNDIR)
	cp examples/RMG/hexadiene/RMG.sh $(RUNDIR)
	export RMG=$(CURDIR); cd $(RUNDIR); ./RMG.sh;
	echo "Results saved to $(RUNDIR)"

# Run all test cases
test_all: 
	export RMG=$(CURDIR); scripts/test_all.sh;

help:
	@echo ""
	@echo "This makefile can be used to build all of the code required by RMG."
	@echo ""
	@echo "Typing \`make\` with no arguments will make all of the executables that come"
	@echo "completely bundled with RMG (indicated in the list below with an asterisk)."
	@echo ""
	@echo "Typing 'make test' will run an example file in the folder 'run'."
	@echo ""
	@echo "Typing \`make clean\` will remove all object files, modules, and executables"
	@echo "for all of the applications, and the 'run' folder."
	@echo ""
	@echo "Otherwise individual dependencies can be specified using \`make <target>\`"
	@echo "where <target> is one of:"
	@echo ""
	@echo "*    fame      to make fame, the pressure dependence code"
	@echo "*    frankie   to make frankie, the density of states code"
	@echo "*    dassl     to make dassl, the basic differential equation solver"
	@echo "     daspk     to make daspk, the diff. eq. solver with sensitivity"
	@echo "*    GATPFit   to make GATPFit"
	@echo "*    RMG       to make RMG"
	@echo "     inchi     to get and make InChI support (attempts to download from iupac)"
	@echo "     symmetry  to get and make SYMMETRY (attempts to download from ucalgary.ca)"
	@echo ""
