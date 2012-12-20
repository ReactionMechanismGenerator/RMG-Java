@echo off

if not defined RMG (
echo Error: RMG variable not defined; please set to the location of your RMG installation.
goto end
)

mkdir "%rmg%"\bin

echo Compiling GATPFit...
cd "%rmg%"\source\GATPFit
g95 -Wall -fbounds-check *.f -L"%rmg%"\bin -lblas -llapack -o "%rmg%"\bin\GATPFit.exe

IF EXIST  "%RMG%"\bin\GATPFit.exe (
        echo GATPFit.exe created succesfully.
) ELSE (
	echo Compilation failed. Make sure blas.dll and lapack.dll are in %RMG%\bin folder. You can download from https://github.com/GreenGroup/RMG-Java/downloads
	GOTO end	
)

echo Compiling fame...
cd "%rmg%"\source\fame
g95 -Wall -fbounds-check -O3 math.f90 states.f90 specfun.f90 _modes.f90 network.f90 io.f90 mastereqn.f90 msc.f90 rs.f90 model.f90 fame.f90 -L"%rmg%"\bin -lblas -llapack -o "%rmg%"\bin\fame.exe

IF EXIST  "%RMG%"\bin\fame.exe (
        echo GATPFit.exe created succesfully.
) ELSE (
	echo Compilation failed. Check g95 output.
	GOTO end	
)

echo Compiling frankie...
cd "%rmg%"\source\frankie
g95 -Wall -fbounds-check *.f90 -o "%rmg%"\bin\frankie.exe


IF EXIST  "%RMG%"\bin\frankie.exe (
        echo GATPFit.exe created succesfully.
) ELSE (
	echo Compilation failed. Check g95 output.
	GOTO end	
)


echo Compiling dassl...
cd "%rmg%"\source\dassl
g95 call_dasslAUTO.f90 daux.f ddassl.f dlinpk.f getflux.f reaction_flux.f res.f res_daepack.f -o "%rmg%"\bin\dasslAUTO.exe

IF EXIST  "%RMG%"\bin\dasslAUTO.exe (
        echo GATPFit.exe created succesfully.
) ELSE (
	echo Compilation failed. Check g95 output.
	GOTO end	
)

echo Compiled Fortran code succesfully. You should have GATPFit, fame, frankie, dasslAUTO executables in "%RMG%"\bin folder.
echo This script does NOT compile the Java code; use the command "ant jar" to do so.

cd "%rmg%"

:end
pause
