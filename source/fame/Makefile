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

# The Fortran 90 compiler to use and flags to use when compiling Fortran 90 
# code
F90=g95
F90FLAGS=

LDFLAGS=-lpthread -lblas -llapack

OBJ0=Species.o Isomer.o Reaction.o DensityOfStates.o Network.o Input.o StrongCollision.o ReservoirState.o RateModel.o Output.o fame.o 

OBJ=$(patsubst %.o,$(BUILDDIR)/%.o,$(OBJ0))

################################################################################

$(BINDIR)/fame.exe: $(OBJ)  
	mkdir -p $(BINDIR)
	$(F90) $(LDFLAGS) $(OBJ) -o $(BINDIR)/fame.exe


$(BUILDDIR)/%.o: %.f90
	mkdir -p $(BUILDDIR)
	$(F90) $(F90FLAGS) -c $< -o $@

clean:
	rm -rf $(BUILDDIR)
	rm -rf $(BINDIR)
