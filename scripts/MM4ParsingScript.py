#uses modified cclib 1.00
import sys
sys.path.append(sys.argv[2])#add $RMG/source to the PYTHONPATH so that import statements below work properly
import logging 
from cclib.parser import MM4
myfile=MM4(sys.argv[1])
myfile.logger.setLevel(logging.ERROR) #cf. http://cclib.sourceforge.net/wiki/index.php/Using_cclib#Additional_information
data=myfile.parse()
print data.natom #print line 1
print [i for i in data.atomnos] #this will put a comma between each element
print data.atomcoords[-1] #print the final optimized atom coordinates
print data.scfenergies[-1]/27.2113845 #print the MM4 heat of formation (cclib gives this in eV, but I have converted here to Hartree); conversion value taken from cclib code; compare CODATA 2006 value 27.21138386(68) eV/Hartree
if(sys.argv[3]=='1'):##sys.argv[3]=1 requests that steric energy also be printed before molarmass
    print data.stericenergy/27.2113845 #print the MM4 final steric energy (cclib gives this in eV, but I have converted here to Hartree); conversion value taken from cclib code; compare CODATA 2006 value 27.21138386(68) eV/Hartree
print data.molmass #print the molecular mass (in amu)
if (data.natom > 1):#for species with more than one atom, display vibrational frequencies and rotational constants
    print [i for i in data.vibfreqs]  
    print data.rotcons[-1] #print the final rotational constants