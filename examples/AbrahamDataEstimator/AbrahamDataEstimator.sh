#! /bin/bash

if ($RMG)
    then
    echo "Environment variable \$RMG is not defined.  Please set to the location of your RMG installation."
    exit
fi

if [ ! -e "$RMG/bin/RMG.jar" ]
then
  echo "$RMG/bin/RMG.jar not found; please compile using ant."; 
  exit
fi

if [ ! -e "Abraham_input.txt" ]
then
  echo "Error: Abraham_input.txt not found. Please create an AbrahamDataEstimator input file before running."; 
  exit
fi

echo "Running AbrahamDataEstimator..."
java -Xmx500m -classpath  $RMG/bin/RMG.jar Abraham 2>&1 | tee RMG.log
echo "AbrahamDataEstimator job completed."
