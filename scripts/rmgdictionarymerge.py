# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="User1"
__date__ ="$Mar 7, 2012 9:21:43 PM$"

# arg[1]: master RMG_Dictionary file
# arg[2]: donor RMG_Dictionary file
# arg[3]: CHEMKIN-MFC species_mapping file for merge
# arg[4]: output merged RMG_Dictionary
if __name__ == "__main__":
	import sys
	comma=","
	tab="\t"
	master=open(sys.argv[1], 'r')
	donor=open(sys.argv[2],'r')
	mfc=open(sys.argv[3],'r')
	merge=open(sys.argv[4], 'w')

	#first, copy the master dictionary file
	for line in master:
	    merge.write(line)
	master.close()

	#merge.write("//following species from donor dictionary:\n")
	merge.write("\n")

	#read in the species_mapping file
	mDict={}
	rDict={}
	rename = False
	line = mfc.readline()
	while not line.startswith('Species_in_donor_mechanism'):#get to the main section with duplicated species
	    if line.startswith('Original_name_in_donor_mechanism'):#read renamed section if it exists
		rename = True
	    elif rename:
		if line == '\n':#stop reading once blank line encountered
		    rename = False
		else:
		    split = line.split()
		    rDict[split[0]]=split[1] #old (donor) name is key, new name is value
	    line=mfc.readline()
	while line != '' and line != '\n':#read the duplicated species section until the end of the file/blank line
	    split = line.split()
	    mDict[split[0]]=split[1] #old (donor) name is key, new name is value
	    line=mfc.readline()
	mfc.close()

	#read in the donor dictionary
	line = donor.readline()
	while line != '':
	    while line == '\n':#get to the next species block
		line = donor.readline()
	    name = line.split()[0]#get the species name
	    writeFlag = True
	    if name in rDict:
		merge.write("\n")
		merge.write(rDict[name]+'\n')
	    elif name in mDict: #if the species is already present, we don't need to add it to the dictionary as it is already there
		writeFlag = False
	    else:#copy the species exactly if it isn't mentioned in the rDict or mDict
		merge.write("\n")
		merge.write(line)#already includes "\n"
	    line = donor.readline()
	    while line !='\n' and line != '':
		if writeFlag:
		    merge.write(line)
		line = donor.readline()
	    line = donor.readline()
	donor.close()
	merge.close()
