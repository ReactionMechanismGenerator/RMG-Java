echo off
if "%rmg%"=="" echo The RMG variable is NOT defined!

mkdir "%rmg%"\bin

echo Compiling GATPFit...
cd "%rmg%"\source\GATPFit
g95 *.f -L"%rmg%"\bin -lblas -llapack -o "%rmg%"\bin\GATPFit.exe

echo Compiling fame...
cd "%rmg%"\source\fame
g95 -O3 math.f90 states.f90 _modes.f90 network.f90 io.f90 mastereqn.f90 msc.f90 rs.f90 model.f90 fame.f90 -L"%rmg%"\bin -lblas -llapack -o "%rmg%"\bin\fame.exe

echo Compiling frankie...
cd "%rmg%"\source\frankie
g95 *.f90 -o "%rmg%"\bin\frankie.exe

echo Compiling dassl...
cd "%rmg%"\source\dassl
g95 call_dasslAUTO.f90 daux.f ddassl.f dlinpk.f getflux.f reaction_flux.f res.f res_daepack.f -o "%rmg%"\bin\dasslAUTO.exe

cd "%rmg%"