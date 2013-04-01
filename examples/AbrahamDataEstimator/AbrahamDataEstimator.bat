@echo off

if not defined RMG (
echo Error: RMG variable not defined; please set to the location of your RMG installation.
goto end
)
if not exist "%RMG%\bin\RMG.jar" (
echo Error: RMG.jar not found; please compile using ant.
goto end
)

if not exist Abraham_input.txt (
echo Error: Abraham_input.txt not found. Please create an AbrahamDataEstimator input file before running.
goto end
)

echo Running AbrahamDataEstimator...
java -Xmx500m -classpath "%RMG%\bin\RMG.jar" Abraham  > RMG.log 2>&1 &
echo AbrahamDataEstimator job completed.

:end
pause
