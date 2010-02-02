cd "%rmg%"\source\GATPFit\src
g95 -o GATPFit *.f
mkdir "%rmg%"\software\GATPFit
copy GATPFit.exe "%rmg%"\software\GATPFit\GATPFit.exe
cd "%rmg%"\source\fame
g95 -o fame src\blas\*.f src\lapack\double\*.f src\lapack\util\*.f src/Simulation.f90 src/Isomer.f90 src/Reaction.f90 src/Input.f90 src/DensityOfStates.f90 src/MasterEqn.f90 src/StrongCollision.f90 src/ReservoirState.f90 src/RateModel.f90 src/fame.f90
mkdir "%rmg%"\software\fame
copy fame.exe "%rmg%"\software\fame\fame.exe
cd "%rmg%"\source\frankie\src
g95 -o frankie *.f90
mkdir "%rmg%"\software\frankie
copy frankie.exe "%rmg%"\software\frankie\frankie.exe
cd "%rmg%"\source\ODESolver\dassl
g95 -o dasslAUTO *.f *.f90
mkdir "%rmg%"\software\ODESolver
copy dasslAUTO.exe "%rmg%"\software\ODESolver\dasslAUTO.exe


