#!/bin/bash

## This base $RMG variable must be set correctly.
# export RMG=/path/to/rmg/

# All the other variables here have default values, as shown:


# Now run RMG on the condition.txt input file.
java -Xmx500m -jar $RMG/bin/RMG.jar condition.txt
