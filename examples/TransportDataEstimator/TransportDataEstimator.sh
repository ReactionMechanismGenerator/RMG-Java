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

if [ ! -e "input.txt" ]
then
  echo "Error: input.txt not found. Please create an TransportDataEstimator input file before running."; 
  exit
fi

echo "Running TransportDataEstimator..."
java -Xmx500m -classpath  $RMG/bin/RMG.jar TransportDataEstimator input.txt 2>&1 | tee RMG.log
echo "TransportDataEstimator job completed."
