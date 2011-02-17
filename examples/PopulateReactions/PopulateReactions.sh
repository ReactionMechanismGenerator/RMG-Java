#! /bin/bash
echo "Running PopulateReactions..."
java -Xmx500m -classpath $RMG/bin/RMG.jar PopulateReactions input.txt 2>&1 | tee RMG.log
echo "PopulateReactions job completed. Results saved in RMG.log"
