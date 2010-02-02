################################################################################
#
#	Makefile for fame
#
################################################################################

# The directory in which the source files can be found
$SOURCEDIR=.

# The directory in which to place temporary compiled files
$BUILDDIR=../../build/frankie

# The directory in which to place compiled executables and JAR files
$BINDIR=../../bin

# The Fortran 90 compiler to use and flags to use when compiling Fortran 90 
# code
F90=g95
F90FLAGS=

LDFLAGS=-lpthread

OBJ=dqed.o calc_freq_code.o

OBJ=$(patsubst %.o,$(BUILDDIR)/%.o,$(OBJ))

################################################################################

$(BINDIR)/frankie.exe: dirs $(OBJ)  
	$(F90) $(LDFLAGS) $(OBJ) -o $(BINDIR)/frankie.exe

dirs:
	mkdir -p $(BUILDDIR)
	mkdir -p $(BINDIR)

$(BUILDDIR)/%.o: %.f90
	$(F90) $(F90FLAGS) -c $< -o $@

clean:
	rm -rf $(BUILDDIR)
	rm -rf $(BINDIR)
