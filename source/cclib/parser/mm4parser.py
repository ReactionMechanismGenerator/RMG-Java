"""
gmagoon 05/03/10: new class for MM4 parsing, based on mopacparser.py, which, in turn, is based on gaussianparser.py from cclib, described below:
cclib (http://cclib.sf.net) is (c) 2006, the cclib development team
and licensed under the LGPL (http://www.gnu.org/copyleft/lgpl.html).
"""

__revision__ = "$Revision: 814 $"


#import re

import numpy

import utils
import logfileparser


def symbol2int(symbol):
    t = utils.PeriodicTable()
    return t.number[symbol]

class MM4(logfileparser.Logfile):
    """An MM4 output file."""

    def __init__(self, *args, **kwargs):

        # Call the __init__ method of the superclass
        super(MM4, self).__init__(logname="MM4", *args, **kwargs)
        
    def __str__(self):
        """Return a string representation of the object."""
        return "MM4 log file %s" % (self.filename)

    def __repr__(self):
        """Return a representation of the object."""
        return 'MM4("%s")' % (self.filename)

    def extract(self, inputfile, line):
        """Extract information from the file object inputfile."""
        
        # Number of atoms.
        # Example:          THE COORDINATES OF    20 ATOMS ARE READ IN.
        if line[0:28] == '          THE COORDINATES OF':

            self.updateprogress(inputfile, "Attributes", self.fupdate)       
            natom = int(line.split()[-5]) #fifth to last component should be number of atoms
            if hasattr(self, "natom"):
                assert self.natom == natom
            else:
                self.natom = natom
        
        # Extract the atomic numbers and coordinates from the optimized (final) geometry
        
        # Example:
#	      FINAL ATOMIC COORDINATE
#           ATOM          X           Y           Z      TYPE
#         C(    1)    -3.21470    -0.22058     0.00000   (  1)
#         H(    2)    -3.30991    -0.87175     0.89724   (  5)
#         H(    3)    -3.30991    -0.87174    -0.89724   (  5)
#         H(    4)    -4.08456     0.47380     0.00000   (  5)
#         C(    5)    -1.88672     0.54893     0.00000   (  1)
#         H(    6)    -1.84759     1.21197    -0.89488   (  5)
#         H(    7)    -1.84759     1.21197     0.89488   (  5)
#         C(    8)    -0.66560    -0.38447     0.00000   (  1)
#         H(    9)    -0.70910    -1.04707    -0.89471   (  5)
#         H(   10)    -0.70910    -1.04707     0.89471   (  5)
#         C(   11)     0.66560     0.38447     0.00000   (  1)
#         H(   12)     0.70910     1.04707     0.89471   (  5)
#         H(   13)     0.70910     1.04707    -0.89471   (  5)
#         C(   14)     1.88672    -0.54893     0.00000   (  1)
#         H(   15)     1.84759    -1.21197    -0.89488   (  5)
#         H(   16)     1.84759    -1.21197     0.89488   (  5)
#         C(   17)     3.21470     0.22058     0.00000   (  1)
#         H(   18)     3.30991     0.87174     0.89724   (  5)
#         H(   19)     4.08456    -0.47380     0.00000   (  5)
#         H(   20)     3.30991     0.87175    -0.89724   (  5)

        if line[0:29] == '      FINAL ATOMIC COORDINATE':


            self.updateprogress(inputfile, "Attributes", self.cupdate)
                    
            self.inputcoords = []
            self.inputatoms = []
            
            headerline = inputfile.next()
            
            atomcoords = []
            line = inputfile.next()
            while len(line.split()) > 0:
                broken = line.split()
                self.inputatoms.append(symbol2int(line[0:10].trim()))
                xc = float(line[17:29])
                yc = float(line[29:41])
                zc = float(line[41:53])
                atomcoords.append([xc,yc,zc])
                line = inputfile.next()

            self.inputcoords.append(atomcoords)

            if not hasattr(self, "natom"):
                self.atomnos = numpy.array(self.inputatoms, 'i')
                self.natom = len(self.atomnos)

#read energy (in kcal/mol, converted to eV)
#       Example:     HEAT OF FORMATION (HFN) AT  298.2 K       =       -42.51 KCAL/MOLE
        if line[0:32] == '     HEAT OF FORMATION (HFN) AT':
            if not hasattr(self, "scfenergies"):
                self.scfenergies = []
            self.scfenergies.append(utils.convertor(self.float(line.split()[-2])/627.5095, "hartree", "eV")) #note conversion from kcal/mol to hartree

        #molecular mass parsing (units will be amu); note that this can occur multiple times in the file, but all values should be the same
        #Example:               FORMULA WEIGHT   :     86.112
        if line[0:33] == '               FORMULA WEIGHT   :':
            self.updateprogress(inputfile, "Attributes", self.fupdate)
	    molmass = self.float(line.split()[-1])
	    if hasattr(self, "molmass"):
                assert self.molmass == molmass #check that subsequent occurences match the original value
            else:
                self.molmass = molmass
        
	  #rotational constants (converted to GHZ)
        #Example:

#          ROTATIONAL CONSTANTS IN CM(-1)
#
#          A =    0.01757641   B =    0.00739763   C =    0.00712013
        #could also read in moment of inertia, but this should just differ by a constant: rot cons= h/(8*Pi^2*I)
        #note that the last occurence of this in the thermochemistry section has reduced precision, so we will want to use the 2nd to last instance
        if line[0:40] == '          ROTATIONAL CONSTANTS IN CM(-1)':
	    blankline = inputfile.next();
            rotinfo=inputfile.next();
            if not hasattr(self, "rotcons"):
                self.rotcons = []
            broken = rotinfo.split()
            sol = 29.9792458 #speed of light in vacuum in 10^9 cm/s, cf. http://physics.nist.gov/cgi-bin/cuu/Value?c|search_for=universal_in!
            a = float(broken[2])*sol 
            b = float(broken[5])*sol
            c = float(broken[8])*sol
            self.rotcons.append([a, b, c]) 

        # Start of the IR/Raman frequency section.
#Example:

        if line[1:10] == 'VIBRATION':
	    line = inputfile.next()
            self.updateprogress(inputfile, "Frequency Information", self.fupdate)
      
            if not hasattr(self, 'vibfreqs'):
                self.vibfreqs = []
            freq = self.float(line.split()[1])
            #self.vibfreqs.extend(freqs)
            self.vibfreqs.append(freq)


if __name__ == "__main__":
    import doctest, mm4parser
    doctest.testmod(mm4parser, verbose=False)
