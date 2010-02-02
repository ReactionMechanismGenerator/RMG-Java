# Default is to build those that come complete with the RMG distribution
base: fame frankie GATPFit dassl

all: fame frankie GATPFit dassl daspk

fame: software/fame/fame.exe

frankie: software/frankie/frankie.exe

dassl: software/ODESolver/dasslAUTO.exe

daspk: software/ODESolver/daspkAUTO.exe

GATPFit: software/GATPFit/GATPFit.exe

software/fame/fame.exe:
	make -C fortran/fame/src/

software/frankie/frankie.exe:
	make -C fortran/frankie/src/

software/ODESolver/dasslAUTO.exe:
	make -C fortran/ODESolver/dassl/

software/ODESolver/daspkAUTO.exe:
	make -C fortran/ODESolver/daspk

software/GATPFit/GATPFit.exe:
	make -C fortran/GATPFit

clean:
	make -C fortran/fame/src/ clean
	make -C fortran/frankie/src/ clean
	make -C fortran/GATPFit clean
	make -C fortran/ODESolver/dassl clean
	make -C fortran/ODESolver/DASPK clean
	rm -f ../software/fame/fame.exe
	rm -f ../software/frankie/frankie.exe
	rm -f ../software/GATPFit/GATPFit.exe
	rm -f ../software/ODESolver/*.exe
