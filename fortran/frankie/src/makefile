# ==============================================================================
#
# 	Makefile for masterEqn
#
# ==============================================================================

# PARAMETERS
# ----------

# Name of compiler
FORTRAN = g95

# Name of executable
BIN = ../frankie.exe

# List of object files in project
OBJ = dqed.o calc_freq_code.o 

# Compiler flags
CFLAGS = -ftrace=full

# Linker flags
LFLAGS = -L/usr/lib -lpthread

# ==============================================================================

# RULES
# -----

# "all" rule
all: $(BIN)

# "clean" rule
clean: 
	rm -rf $(OBJ) $(BIN) *.mod

# Linker rule
$(BIN): $(OBJ)
	$(FORTRAN) $(OBJ) -o $(BIN) $(LFLAGS)

# Compiler rule
%.o: %.f90
	$(FORTRAN) $(CFLAGS) -c $< -o $@
