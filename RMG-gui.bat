@echo off

if not defined RMG (
echo Error: RMG environment variable not defined; please set to the location of your RMG installation.
goto end
)
if not exist %RMG%\bin\RMG.jar (
echo Error: RMG.jar not found; please compile using ant.
goto end
)

echo Starting RMG GUI...
java -Xmx500m -classpath %RMG%\bin\RMG.jar GUI

:end
pause
