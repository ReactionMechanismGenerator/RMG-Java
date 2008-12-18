# ==============================================================================
#
# 	Makefile for masterEqn
#
# ==============================================================================

# PARAMETERS
# ----------

# Name of compiler
FORTRAN = g95
#FORTRAN = gfortran

# Name of executable
BIN = ../fame.exe

# List of object files in project
OBJ = fileIO.o Simulation.o Species.o Reaction.o DensityOfStates.o FastEGME.o FullEGME.o RateModel.o fame.o 

# Compiler flags
CFLAGS = -ftrace=full
#CFLAGS =

# Linker flags
LFLAGS = -L/usr/lib -lpthread -lblas -llapack

# ==============================================================================

# RULES
# -----

# "all" rule
all: $(BIN)

# "clean" rule
clean: 
	rm -rf $(OBJ) $(BIN) *.mod ../fame_output.txt

# Linker rule
$(BIN): $(OBJ)
	$(FORTRAN) $(OBJ) -o $(BIN) $(LFLAGS)
	touch ../fame_output.txt

# Compiler rule
%.o: %.f90
	$(FORTRAN) $(CFLAGS) -c $< -o $@
