#!/bin/bash

## This base $RMG variable must be set correctly.
# export RMG=/path/to/rmg/

# All the other variables here have default values, as shown:

# Location of databases (which one to use is defined in the condition file)
export RMG_DATABASES=$RMG/databases

# Scratch directory, unique to this job (be careful on shared computers).
# Fast I/O is beneficial.
# Default is the current working directory from which you run RMG.
export RMG_JOB_SCRATCH=$PWD

# Output directory, unique to this job (be careful on shared computers).
# This is where output you want to keep are saved.
# Default is the current working directory from which you run RMG.
export RMG_JOB_OUTPUT=$PWD

# Now run RMG on the condition.txt input file.
java -Xmx500m -jar $RMG/bin/RMG.jar condition.txt
