#uses MoleCoor
#argv[1] = input file path
#argv[2] = output file name
#argv[3] = molecule name
#argv[4] = $RMG/source/MoleCoor path
import sys
sys.path.append(sys.argv[4])#add $RMG/source/MoleCoor to the PYTHONPATH so that import statements below work properly
import MolecularCoordinates
mg = MolecularCoordinates.readMOLFileWithConnectivity(sys.argv[1])
mg.writeMM4File(sys.argv[2],sys.argv[3])