################################################################################
#
#	Makefile for fame
#
################################################################################

# The directory in which the source files can be found
SOURCEDIR=.

# The directory in which to place temporary compiled files
BUILDDIR=../../build/fame

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
  F90_FAME_FLAGS = -r8 # used for FAME
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

  F90_FAME_FLAGS = -fdefault-real-8 # used for FAME
  F90FLAGS_NDEBUG = $(F90FLAGS) # used for dassl and daspk
endif ###### END OF gfortran SETTINGS

LDFLAGS=$(F90_EXTRA_LDFLAGS)

LIBS=-lblas -llapack
# if you built ATLAS you'll want this instead: 
# LIBS=-latlas -lptf77blas -lptcblas -llapack

OBJ0=math.o specfun.o states.o _modes.o network.o io.o mastereqn.o msc.o rs.o model.o fame.o

OBJ=$(patsubst %.o,$(BUILDDIR)/%.o,$(OBJ0))

################################################################################

$(BINDIR)/fame.exe: $(OBJ)  
	mkdir -p $(BINDIR)
	$(F90) $(F90_FAME_FLAGS) $(LDFLAGS) -o $(BINDIR)/fame.exe $(OBJ) $(LIBS)

$(BUILDDIR)/%.o: %.f90
	mkdir -p $(BUILDDIR)
	$(F90) $(F90_FAME_FLAGS) $(F90FLAGS) -c $< -o $@

clean:
	rm -rf $(BUILDDIR)
	rm -rf $(BINDIR)/fame.exe
