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
OBJ = Species.o Isomer.o Reaction.o DensityOfStates.o Network.o Input.o StrongCollision.o ReservoirState.o RateModel.o Output.o fame.o 

# Compiler flags
CFLAGS = -ftrace=full -fbounds-check -Wall

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
#	touch ../fame_output.txt
	cp $(BIN) $(RMG)/software/fame/fame.exe

# Compiler rule
%.o: %.f90
	$(FORTRAN) $(CFLAGS) -c $< -o $@
