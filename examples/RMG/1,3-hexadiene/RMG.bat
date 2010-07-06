@echo off

REM RMG 3.2 Windows batch script for RMG execution
REM Put me in the directory containing the condition file and double-click to run RMG.
REM This assumes that the condition file is called condition.txt.
REM This script will check that the necessary setup has been done before running RMG.
REM Output from RMG will be logged to the file RMG.log.

if not defined RMG (
echo Error: RMG environment variable not defined; please set to the location of your RMG installation.
goto end
)
if not exist %RMG%\bin\RMG.jar (
echo Error: RMG.jar not found; please compile using ant.
goto end
)
if not exist %RMG%\bin\blas.dll (
echo Error: blas.dll not found. A 32-bit version can be downloaded from http://github.com/GreenGroup/RMG-Java/downloads.
goto end
)
if not exist %RMG%\bin\lapack.dll (
echo Error: lapack.dll not found. A 32-bit version can be downloaded from http://github.com/GreenGroup/RMG-Java/downloads.
goto end
)
if not exist %RMG%\bin\dasslAUTO.exe (
echo Error: dasslAUTO.exe not found. Please run make.bat to compile DASSL.
goto end
)
if not exist %RMG%\bin\GATPFit.exe (
echo Error: GATPFit.exe not found. Please run make.bat to compile GATPFit.
goto end
)
if not exist %RMG%\bin\fame.exe (
echo Error: fame.exe not found. Please run make.bat to compile fame.
goto end
)
if not exist %RMG%\bin\frankie.exe (
echo Error: frankie.exe not found. Please run make.bat to compile frankie.
goto end
)

if not exist condition.txt (
echo Error: condition.txt not found. Please create an RMG condition file before running.
goto end
)

echo Running RMG...
java -Xmx500m -classpath %RMG%\bin\RMG.jar RMG condition.txt > RMG.log 2>&1 &
echo RMG job completed.

:end
pause
