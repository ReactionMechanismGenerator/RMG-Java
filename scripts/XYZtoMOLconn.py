#uses MoleCoor
#argv[1] = input (XYZ) file path
#argv[2] = output (MOL) file name
#argv[3] = molecule name
#argv[4] = tolerance for connectivity perception
#argv[5] = $RMG/source/MoleCoor path
import sys
sys.path.insert(1, sys.argv[5])#add $RMG/source/MoleCoor to the PYTHONPATH so that import statements below work properly
import MolecularCoordinates
mg = MolecularCoordinates.readXYZFile(sys.argv[1])
mg.perceiveConnectivity(tol=float(sys.argv[4]))
mg.writeMOLFile(sys.argv[2],sys.argv[3], connectivity=True)
