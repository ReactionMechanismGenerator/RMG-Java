#!/bin/bash
if [ "x$RMG" = "x" ]; then
   echo "You should set RMG variable to the main RMG directory"
   exit
fi
echo RMG variable is set to $RMG
RUNDIR=$RMG/run
mkdir -p $RUNDIR
REPORT=$RUNDIR/report.rmg
echo This script will run the examples in $RMG/examples/RMG directory  
echo It takes about ten minutes and outputs are stored in $RUNDIR
echo RMG test calculations started on `date`      | tee $REPORT
for INPUT in {hexadiene,cyclopropane_QM,minimal,liquidphase,butane_pruning,diethylsulfide}
  do
     cd $RUNDIR
     mkdir -p $INPUT
     cp $RMG/examples/RMG/$INPUT/condition.txt $INPUT
     cd $INPUT
     echo "*********" Running $INPUT  "********"  | tee $REPORT
     java -jar $RMG/bin/RMG.jar condition.txt > out.rmg
     OUTPUT=$RUNDIR/$INPUT/RMG.log
     if grep -q "Model generation completed" $OUTPUT ; then
        echo "Model generation completed"         | tee $REPORT
     else
        echo "Model generation failed"            | tee $REPORT
     fi 
     grep "The model core has" $OUTPUT | tail -1  | tee $REPORT
     grep "The model edge has" $OUTPUT | tail -1  | tee $REPORT
     grep "Running time"       $OUTPUT | tail -1  | tee $REPORT
     grep "Memory used"        $OUTPUT | tail -1  | tee $REPORT
     grep "ERROR:"            $OUTPUT             | tee $REPORT
     grep "CRITICAL"          $OUTPUT             | tee $REPORT
  done
echo RMG test calculations finished on `date`     | tee $REPORT
