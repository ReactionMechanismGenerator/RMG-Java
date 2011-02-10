@echo off

REM RMG 3.2 Windows batch script for CheckForwardAndReverseRateCoefficient execution
REM Put me in the directory containing the chem.inp file and double-click to run.
REM This assumes that the input file is called chem.inp.
REM This script will check that the necessary setup has been done before running class.
REM Output from class will be logged to the file RMG.log.

if not defined RMG (
echo Error: RMG environment variable not defined; please set to the location of your RMG installation.
goto end
)
if not exist "%RMG%\bin\RMG.jar" (
echo Error: RMG.jar not found; please compile using ant.
goto end
)

if not exist chem.inp (
echo Error: chem.inp not found. Please create an input file before running.
goto end
)

echo Running CheckForwardAndReverseRateCoefficients...
java -Xmx500m -classpath "%RMG%\bin\RMG.jar" CheckForwardAndReverseRateCoefficients chem.inp > RMG.log 2>&1 &
echo Job completed.

:end
pause
