################################################################################
#
#	Makefile for frankie
#
################################################################################
.PHONY: clean

# The directory in which the source files can be found
SOURCEDIR=.

# The directory in which to place temporary compiled files
BUILDDIR=../../build/frankie

# The directory in which to place compiled executables and JAR files
BINDIR=../../bin

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
  # if gfortran>4.3 then add -fbacktrace (it's not supported in earlier versions)
  ifeq ($(shell expr "$(f90version)" : '[4-9]\.[3-9]'),3)
    F90FLAGS += -fbacktrace
  endif
  
  # detect buggy versions of gfortran
  # See http://gcc.gnu.org/bugzilla/show_bug.cgi?id=47694
  BADGFORTRAN =
  ifeq ($(shell expr "$(gfortranversion)" : '4\.[3-5]'),3) # 3.3 - 4.5, potentially buggy
    BADGFORTRAN = OhNoItsABadOne
    ifeq ($(shell expr "$(gfortranversion)" : '4\.4\.[6-9]'),5) # 4.4.6+ OK
      BADGFORTRAN =
    endif
    ifeq ($(shell expr "$(gfortranversion)" : '4\.5\.[3-9]'),5) # 4.5.3+ OK
      BADGFORTRAN =
    endif
  endif
  ifdef BADGFORTRAN
    GFORTRANTEST = echo "You're using gfortran version $(gfortranversion) which contains\
    a bug that breaks this program. Please change to a different version.\
    One of 4.4.6+, 4.5.3+, 4.6.0+, or before 4.3.0 should work. See\
    https://github.com/GreenGroup/RMG-Java/issues/257 for details."; false
  endif
  
  # Call as `make MACOS=true F90=gfortran` if you want MacOS X 10.6+ 64-bit intel core2 features
  ifdef MACOS
    F90FLAGS += -arch x86_64 -march=core2
    F90_EXTRA_LDFLAGS +=  -framework vecLIB
  endif
  F90FLAGS_NDEBUG = $(F90FLAGS) # used for dassl and daspk
endif ###### END OF gfortran SETTINGS

LDFLAGS= $(F90_EXTRA_LDFLAGS)

LIBS=

F90_EXTRA_FLAGS = #in case the user needs to add anything, eg...
#F90_EXTRA_FLAGS='-I""$(BUILDDIR)""' # for early versions of gfortran

OBJ0=dqed.o calc_freq_code.o

OBJ=$(patsubst %.o,$(BUILDDIR)/%.o,$(OBJ0))

################################################################################

$(BINDIR)/frankie.exe: $(OBJ)
	$(GFORTRANTEST)
	mkdir -p $(BINDIR)
	$(F90) $(LDFLAGS) -o $(BINDIR)/frankie.exe $(OBJ) $(LIBS)

$(BUILDDIR)/%.o: %.f90
	$(GFORTRANTEST)
	mkdir -p $(BUILDDIR)
	$(F90) $(F90FLAGS) $(F90_EXTRA_FLAGS) -c $< -o $@

clean:
	rm -rf $(BUILDDIR)
	rm -rf $(BINDIR)/frankie.exe
