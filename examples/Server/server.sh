#! /bin/bash
echo "Running PopulateReactionsServer..."
java -Xmx500m -classpath $RMG/bin/RMG.jar PopulateReactionsServer input.txt 2>&1
#java -Xmx500m -classpath $RMG/build/RMG/ PopulateReactionsServer input.txt 2>&1
echo "PopulateReactions job completed. Results saved in RMG.log"
