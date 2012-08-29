# gmagoon 7/22/09: writes the lowest energy conformation for mole file in
# arg#1 to mole file in arg#2 (and corresponding crude mole file in arg#3), based on UFF energy of arg#4 embeddings
# (optimized); should be reproducible due to use of consistent randomSeed
#updated 8/12/09 for Q2 2009 version of RDKit (rdkit "import" lines of script)
#arg#5 should contain the absolute path of the RDBASE environment variable
import sys
sys.path.insert(1, sys.argv[5])#add $RDBASE to the PYTHONPATH so that import statements below work properly
from rdkit import Chem
from rdkit.Chem import AllChem
attempts=int(sys.argv[4])
m = Chem.MolFromMolFile(sys.argv[1], removeHs=False)
#m2=Chem.AddHs(m)
AllChem.EmbedMultipleConfs(m, attempts,randomSeed=1)
m2crude = Chem.Mol(m.ToBinary()) #make a copy of the (crude) coordinates via ToBinary
energy=0.0
minEid=0;
lowestE=9.999999e99;#start with a very high number, which would never be reached
for i in range(m.GetNumConformers()):
    AllChem.UFFOptimizeMolecule(m,confId=i)
    energy=AllChem.UFFGetMoleculeForceField(m,confId=i).CalcEnergy()
    if (energy < lowestE):
        minEid = i
        lowestE = energy
    #energy.append(AllChem.UFFGetMoleculeForceField(m,confId=i).CalcEnergy())
f=open(sys.argv[2], 'w')
print >>f,Chem.MolToMolBlock(m,confId=minEid)
f.close()
f=open(sys.argv[3], 'w')#write crude coordinates
print >>f,Chem.MolToMolBlock(m2crude,confId=minEid)
f.close()