#!/bin/bash

## This base $RMG variable must be set correctly.
# export RMG=/path/to/rmg/

# All the other variables here have default values, as shown:

# Location of databases (which one to use is defined in the condition file)
export RMG_DATABASES=$RMG/databases

# Scratch directory, unique to this job (be careful on shared computers).
# Fast I/O is beneficial. 
# Stuff stored here may be useful for debugging but in general you don't need to keep it.
# One option which should be safe on shared computers is:
# export RMG_JOB_SCRATCH=`mktemp -d -t RMG.XXXXXX`
# Or if you are using a PBS or LAVA queuing system you could use the job identifier.
# Default is the current working directory from which you run RMG.
export RMG_JOB_SCRATCH=$PWD

# Output directory, unique to this job (be careful on shared computers).
# This is where output that you want to keep is saved.
# Default is the current working directory from which you run RMG.
export RMG_JOB_OUTPUT=$PWD

# Quantum Mechanics Thermodynamics library directory.
# Default is 'QMThermoLibrary' inside the current job's output directory
# Making it common to all your RMG jobs will allow them to pool all their QM results in the same place
export RMG_QM_LIBRARY=$RMG_JOB_OUTPUT/QMThermoLibrary

# Quantum Mechanics Thermodynamics calculations results directory.
# Default is 'QMfiles' inside the current job's output directory.
# Making it common to all your RMG jobs will allow them to pool all their QM results in the same place
# If you set "KeepQMFiles: no" in your condition file, this directory will be deleted.
export RMG_QM_CALCS=$RMG_JOB_OUTPUT/QMfiles

# Now run RMG on the condition.txt input file.
java -Xmx500m -jar $RMG/bin/RMG.jar condition.txt
