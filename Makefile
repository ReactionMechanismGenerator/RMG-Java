################################################################################
#
#	Makefile for RMG
#
################################################################################

# The directory in which the source files can be found
$SOURCEDIR = source

# The directory in which to place temporary compiled files
$BUILDDIR = build

# The directory in which to place compiled executables and JAR files
$BINDIR = bin

# The Fortran 90 compiler to use and flags to use when compiling Fortran 90 
# code
F90 = g95
F90FLAGS = -fbounds-check -Wall

################################################################################

# Default is to build those that come complete with the RMG distribution
base: dirs fame frankie GATPFit dassl RMG

# You can also build everything
all: dirs fame frankie GATPFit dassl daspk RMG

RMG: dirs $(BINDIR)/RMG.jar

fame: dirs $(BINDIR)/fame.exe

frankie: dirs $(BINDIR)/frankie.exe

dassl: dirs $(BINDIR)/dasslAUTO.exe

daspk: dirs $(BINDIR)/daspkAUTO.exe

GATPFit: dirs $(BINDIR)/GATPFit.exe

$(BINDIR)/RMG.jar:
	ant compile
	ant jar

$(BINDIR)/fame.exe:
	make -C $(SOURCEDIR)/fame

$(BINDIR)/frankie.exe:
	make -C $(SOURCEDIR)/frankie

$(BINDIR)/dasslAUTO.exe:
	make -C $(SOURCEDIR)/dassl

$(BINDIR)/daspkAUTO.exe:
	make -C $(SOURCEDIR)/daspk

$(BINDIR)/GATPFit.exe:
	make -C $(SOURCEDIR)/GATPFit

dirs:
	mkdir -p $(BUILDDIR)
	mkdir -p $(BINDIR)

clean:
	make -C $(SOURCEDIR)/fame clean
	make -C $(SOURCEDIR)/frankie clean
	make -C $(SOURCEDIR)/dassl clean
	make -C $(SOURCEDIR)/daspk clean
	make -C $(SOURCEDIR)/GATPFit clean
	ant clean
	rm -rf $(BUILDDIR)
	rm -rf $(BINDIR)

help:
	@echo ""
	@echo "This makefile can be used to build all of the code required by RMG."
	@echo ""
	@echo "Typing \`make\` with no arguments will make all of the executables that come"
	@echo "completely bundled with RMG (indicated in the list below with an asterisk)."
	@echo ""
	@echo "Typing \`make clean\` will remove all object files, modules, and executables"
	@echo "for all of the applications."
	@echo ""
	@echo "Otherwise individual dependencies can be specified using \`make <target>\`"
	@echo "where <target> is one of:"
	@echo ""
	@echo "*    fame      to make fame, the pressure dependence code"
	@echo "*    frankie   to make frankie, the density of states code"
	@echo "*    dassl     to make dassl, the basic differential equation solver"
	@echo "     daspk     to make daspk, the diff. eq. solver with sensitivity"
	@echo "*    GATPFit   to make GATPFit"
	@echo "*    RMG       to make RMG"
	@echo ""
