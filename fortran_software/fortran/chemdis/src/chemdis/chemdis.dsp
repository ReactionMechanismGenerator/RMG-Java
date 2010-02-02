# Microsoft Developer Studio Project File - Name="chemdis" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Application" 0x0101

CFG=chemdis - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "chemdis.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "chemdis.mak" CFG="chemdis - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "chemdis - Win32 Release" (based on "Win32 (x86) Application")
!MESSAGE "chemdis - Win32 Debug" (based on "Win32 (x86) Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
MTL=midl.exe
RSC=rc.exe

!IF  "$(CFG)" == "chemdis - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE F90 /compile_only /nologo /warn:nofileopt /winapp
# ADD F90 /compile_only /nologo /warn:nofileopt /winapp
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /YX /FD /c
# ADD BASE MTL /nologo /D "NDEBUG" /mktyplib203 /win32
# ADD MTL /nologo /D "NDEBUG" /mktyplib203 /win32
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 version.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:windows /machine:I386
# ADD LINK32 version.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:windows /machine:I386

!ELSEIF  "$(CFG)" == "chemdis - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE F90 /check:bounds /compile_only /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt /winapp
# ADD F90 /check:bounds /compile_only /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt /winapp
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_MBCS" /YX /FD /GZ  /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_MBCS" /YX /FD /GZ  /c
# ADD BASE MTL /nologo /D "_DEBUG" /mktyplib203 /win32
# ADD MTL /nologo /D "_DEBUG" /mktyplib203 /win32
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 version.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:windows /incremental:no /debug /machine:I386 /pdbtype:sept
# ADD LINK32 version.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:windows /incremental:no /debug /machine:I386 /pdbtype:sept

!ENDIF 

# Begin Target

# Name "chemdis - Win32 Release"
# Name "chemdis - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=..\cdaux.f
DEP_F90_CDAUX=\
	"..\cdcontrl.fh"\
	"..\cdisprop.fh"\
	"..\cdnos.fh"\
	"..\cdparams.fh"\
	"..\cdrange.fh"\
	"..\cdrates.fh"\
	"..\cdwell0.fh"\
	"..\cdwell1.fh"\
	"..\cdwell2.fh"\
	"..\straightBS.fh"\
	
# End Source File
# Begin Source File

SOURCE=..\cdbits.f
DEP_F90_CDBIT=\
	"..\cdparams.fh"\
	
# End Source File
# Begin Source File

SOURCE=..\cdchem.f
DEP_F90_CDCHE=\
	"..\cdcolls.fh"\
	"..\cdcontrl.fh"\
	"..\cdisprop.fh"\
	"..\cdparams.fh"\
	"..\cdrange.fh"\
	"..\cdrates.fh"\
	"..\cdwell0.fh"\
	"..\cdwell1.fh"\
	"..\cdwell2.fh"\
	
# End Source File
# Begin Source File

SOURCE=..\cddiss.f
DEP_F90_CDDIS=\
	"..\cdASATempPres.fh"\
	"..\cdcolls.fh"\
	"..\cdcontrl.fh"\
	"..\cdisprop.fh"\
	"..\cdlabels.fh"\
	"..\cdparams.fh"\
	"..\cdrange.fh"\
	"..\cdrates.fh"\
	"..\cdwell0.fh"\
	"..\cdwell1.fh"\
	"..\cdwell2.fh"\
	
# End Source File
# Begin Source File

SOURCE=..\cdffits.f
DEP_F90_CDFFI=\
	"..\cdffit.fh"\
	"..\cdffunc.fh"\
	"..\cdlabels.fh"\
	"..\cdparams.fh"\
	"..\cdrange.fh"\
	"..\cdwell0.fh"\
	
# End Source File
# Begin Source File

SOURCE=..\cdffuns.f
# End Source File
# Begin Source File

SOURCE=..\cdfouts.f
DEP_F90_CDFOU=\
	"..\cdffunc.fh"\
	"..\cdlabels.fh"\
	"..\cdlimfit.fh"\
	"..\cdparams.fh"\
	"..\cdrange.fh"\
	"..\cdrates.fh"\
	"..\cdwell0.fh"\
	
# End Source File
# Begin Source File

SOURCE=..\cdgets.f
DEP_F90_CDGET=\
	"..\cdcheb.fh"\
	"..\cdcolls.fh"\
	"..\cdcontrl.fh"\
	"..\cdffit.fh"\
	"..\cdgamma.fh"\
	"..\cdisprop.fh"\
	"..\cdlabels.fh"\
	"..\cdnos.fh"\
	"..\cdparams.fh"\
	"..\cdrange.fh"\
	"..\cdwell0.fh"\
	"..\cdwell1.fh"\
	"..\cdwell2.fh"\
	"..\straightBS.fh"\
	
# End Source File
# Begin Source File

SOURCE=..\cdkfits.f
DEP_F90_CDKFI=\
	"..\cdlabels.fh"\
	"..\cdlimfit.fh"\
	"..\cdparams.fh"\
	"..\cdrange.fh"\
	"..\cdrates.fh"\
	"..\cdwell0.fh"\
	
# End Source File
# Begin Source File

SOURCE=..\cdkouts.f
DEP_F90_CDKOU=\
	"..\cdcheb.fh"\
	"..\cdcolls.fh"\
	"..\cdcontrl.fh"\
	"..\cdkfit.fh"\
	"..\cdlabels.fh"\
	"..\cdparams.fh"\
	"..\cdrange.fh"\
	"..\cdrates.fh"\
	"..\cdwell0.fh"\
	
# End Source File
# Begin Source File

SOURCE=..\cdrecur.f
DEP_F90_CDREC=\
	"..\cdisprop.fh"\
	"..\cdparams.fh"\
	
# End Source File
# Begin Source File

SOURCE=..\cdsets.f
DEP_F90_CDSET=\
	"..\cdcontrl.fh"\
	"..\cdisprop.fh"\
	"..\cdlabels.fh"\
	"..\cdparams.fh"\
	"..\cdwell0.fh"\
	"..\cdwell1.fh"\
	"..\cdwell2.fh"\
	
# End Source File
# Begin Source File

SOURCE=..\cdstab.f
DEP_F90_CDSTA=\
	"..\cdcolls.fh"\
	"..\cdcontrl.fh"\
	"..\cdisprop.fh"\
	"..\cdparams.fh"\
	
# End Source File
# Begin Source File

SOURCE=..\cdstr.f
# End Source File
# Begin Source File

SOURCE=..\cdsubs.f
DEP_F90_CDSUB=\
	"..\cdgamma.fh"\
	"..\cdisprop.fh"\
	"..\cdnos.fh"\
	"..\cdparams.fh"\
	"..\cdrhovar.fh"\
	"..\cdwell0.fh"\
	"..\straightBS.fh"\
	
# End Source File
# Begin Source File

SOURCE=..\chebyshev.f
DEP_F90_CHEBY=\
	"..\cdcheb.fh"\
	"..\cdparams.fh"\
	"..\cdrange.fh"\
	"..\cdrates.fh"\
	
# End Source File
# Begin Source File

SOURCE=..\chemdis2.f
DEP_F90_CHEMD=\
	"..\cdcheb.fh"\
	"..\cdcontrl.fh"\
	"..\cdgamma.fh"\
	"..\cdlabels.fh"\
	"..\cdlimfit.fh"\
	"..\cdparams.fh"\
	"..\cdrange.fh"\
	"..\cdrates.fh"\
	"..\cdwell0.fh"\
	
# End Source File
# Begin Source File

SOURCE=..\gammln.f
# End Source File
# Begin Source File

SOURCE=..\genRhos.f
DEP_F90_GENRH=\
	"..\cdcontrl.fh"\
	"..\cdisprop.fh"\
	"..\cdnos.fh"\
	"..\cdparams.fh"\
	"..\cdwell0.fh"\
	"..\cdwell1.fh"\
	"..\straightBS.fh"\
	
# End Source File
# Begin Source File

SOURCE=..\PrintBarkerFiles.f
DEP_F90_PRINT=\
	"..\cdcontrl.fh"\
	"..\cdisprop.fh"\
	"..\cdlabels.fh"\
	"..\cdnos.fh"\
	"..\cdparams.fh"\
	"..\cdrange.fh"\
	"..\cdwell0.fh"\
	"..\cdwell1.fh"\
	"..\cdwell2.fh"\
	"..\straightBS.fh"\
	
# End Source File
# Begin Source File

SOURCE=..\PrintSteadyInput.f
DEP_F90_PRINTS=\
	"..\cdcontrl.fh"\
	"..\cdlabels.fh"\
	"..\cdparams.fh"\
	"..\cdrange.fh"\
	"..\cdrates.fh"\
	"..\cdwell0.fh"\
	
# End Source File
# Begin Source File

SOURCE=..\straightBS.f
DEP_F90_STRAI=\
	"..\straightBS.fh"\
	
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl;fi;fd"
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# End Target
# End Project
