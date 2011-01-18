################################################################################
#
#	Makefile for frankie
#
################################################################################

# The directory in which the source files can be found
SOURCEDIR=.

# The directory in which to place temporary compiled files
BUILDDIR=../../build/frankie

# The directory in which to place compiled executables and JAR files
BINDIR=../../bin

# The Fortran 90 compiler to use and flags to use when compiling Fortran 90 
# code
#F90=g95
#F90FLAGS=-fmod=$(BUILDDIR) -ftrace=full

F90=gfortran
# --fpe-trap=underflow,overflow causes problems at least with fame. Maybe a bug in gcc? (Maybe a bug in fame.)
F90FLAGS = -ffpe-trap=invalid,zero -ftrapv -fbounds-check -frange-check \
           -ggdb -J""$(BUILDDIR)"" -O3  -Wall -Wno-unused $(backtrace)
# if gfortran>4.3 then add -fbacktrace (it's not supported in earlier versions)
backtrace = $(shell gfortran --version 2>/dev/zero|grep -iqs '^GNU Fortran.* [4-9]\.[3-9]\.[0-9]' && echo "-fbacktrace")

LDFLAGS= $(F90_EXTRA_LDFLAGS) -lpthread

OBJ0=dqed.o calc_freq_code.o

OBJ=$(patsubst %.o,$(BUILDDIR)/%.o,$(OBJ0))

################################################################################

$(BINDIR)/frankie.exe: $(OBJ)  
	mkdir -p $(BINDIR)
	$(F90) $(LDFLAGS) $(OBJ) -o $(BINDIR)/frankie.exe

$(BUILDDIR)/%.o: %.f90
	mkdir -p $(BUILDDIR)
	$(F90) $(F90FLAGS) -c $< -o $@

clean:
	rm -rf $(BUILDDIR)
	rm -rf $(BINDIR)/frankie.exe
