if "%rmg%"=="" echo The RMG variable is NOT defined!

mkdir "%rmg%"\bin

echo Compiling GATPFit...
cd "%rmg%"\source\GATPFit
g95 *.f blas\*.f lapack\double\*.f lapack\util\*.f -o "%rmg%"\bin\GATPFit.exe

echo Compiling fame...
cd "%rmg%"\source\fame
g95 Species.f90 Isomer.f90 Reaction.f90 DensityOfStates.f90 Network.f90 Input.f90 StrongCollision.f90 ReservoirState.f90 RateModel.f90 Output.f90 fame.f90 blas\*.f lapack\double\*.f lapack\util\*.f -o "%rmg%"\bin\fame.exe

echo Compiling frankie...
cd "%rmg%"\source\frankie
g95 *.f90 -o "%rmg%"\bin\frankie.exe

echo Compiling dassl...
cd "%rmg%"\source\dassl
g95 *.f *.f90 -o "%rmg%"\bin\dasslAUTO.exe

