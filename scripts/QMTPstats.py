# arg[1]: RMG log file from a run with QM with the "both" option
if __name__ == "__main__":
	import sys

#*****Attempt #5 on species VLEXMQNIVAQEHV-UHFFFAOYAI (InChI=1/C10H14/c1-2-3-4-7-10-8-5-6-9-10/h2,7,10H,5-6,8-9H2,1H3) failed. Will attempt a new keyword.
#*****Imaginary freqencies found:
#*****Attempt #6 on species VLEXMQNIVAQEHV-UHFFFAOYAI (InChI=1/C10H14/c1-2-3-4-7-10-8-5-6-9-10/h2,7,10H,5-6,8-9H2,1H3) failed. Will attempt a new keyword.
#*****Imaginary freqencies found:
#*****Attempt #7 on species VLEXMQNIVAQEHV-UHFFFAOYAI (InChI=1/C10H14/c1-2-3-4-7-10-8-5-6-9-10/h2,7,10H,5-6,8-9H2,1H3) failed. Will attempt a new keyword.
#*****Imaginary freqencies found:
#*****Attempt #8 on species VLEXMQNIVAQEHV-UHFFFAOYAI (InChI=1/C10H14/c1-2-3-4-7-10-8-5-6-9-10/h2,7,10H,5-6,8-9H2,1H3) failed. Will attempt a new keyword.
#*****Imaginary freqencies found:
#*****Attempt #9 on species VLEXMQNIVAQEHV-UHFFFAOYAI (InChI=1/C10H14/c1-2-3-4-7-10-8-5-6-9-10/h2,7,10H,5-6,8-9H2,1H3) failed. Will attempt a new keyword.
#*****Imaginary freqencies found:
#*****Final MOPAC attempt (#10) on species VLEXMQNIVAQEHV-UHFFFAOYAI (InChI=1/C10H14/c1-2-3-4-7-10-8-5-6-9-10/h2,7,10H,5-6,8-9H2,1H3) failed. Trying to use Gaussian.
#*****Attempt #0 on species VLEXMQNIVAQEHV-UHFFFAOYAI (InChI=1/C10H14/c1-2-3-4-7-10-8-5-6-9-10/h2,7,10H,5-6,8-9H2,1H3) failed. Will attempt a new keyword.
#Attempt #1 on species VLEXMQNIVAQEHV-UHFFFAOYAI (InChI=1/C10H14/c1-2-3-4-7-10-8-5-6-9-10/h2,7,10H,5-6,8-9H2,1H3) succeeded.
#Point group: Cs
#Thermo for VLEXMQNIVAQEHV-UHFFFAOYAI: 48.57     112.68  41.61   54.59   66.16   75.82   90.6    101.25  117.34
#HBI-based thermo for OSUAYNIVCVOFCQ-UHFFFAOYATmult3(InChI=1/C10H12/c1-2-3-4-7-10-8-5-6-9-10/h2,7H,1,5-6,8-9H2/mult3): 115.97    106.43      39.2    51.65   62.64   71.73   85.48   95.24   110.3
#Attempt #1 on species DVGVDWBLASWERI-UHFFFAOYAM (InChI=1/C10H14/c1-2-3-4-7-10-8-5-6-9-10/h2,10H,1,5-9H2) succeeded.
#Point group: C1
#Thermo for DVGVDWBLASWERI-UHFFFAOYAM: 39.96     104.61  40.87   54.22   66.0    75.78   90.65   101.32  117.41
#HBI-based thermo for OSUAYNIVCVOFCQ-UHFFFAOYATmult3(InChI=1/C10H12/c1-2-3-4-7-10-8-5-6-9-10/h2,7H,1,5-6,8-9H2/mult3): 119.26    109.4       39.5    50.54   60.7    69.44   82.99   92.81   108.53
	MopacKeywords = 10
	GaussianKeywords = 36

	#initialize counters

	preexistMopac=0
	preexistGaussian=0
	preexistMopacList=[]
	preexistGaussianList=[]
	masterList=[]
	connMismatchList=[]
	newMopacArray=[0 for x in range(MopacKeywords)]
	newGaussianArray=[0 for x in range(GaussianKeywords)]
	newFail=0
	connMismatchSpeciesRetry=0
	connMismatchFail=0
	#read in the file
	iin=open(sys.argv[1], 'r')
	line = iin.readline()
	while line != '':#'' marks end of file
	    if ('Attempt #1' in line):#the beginning of a block for a new attempt
		#print line
		connMismatchFlag=False
		attemptCounter=1
		while (not line.startswith('Attempt') and not line.startswith('*****Final attempt')):#this loop gets to the successful attempt or the last attempt, counting the number of intermediate failed attempts in the process
		    line = iin.readline()
		    if line.startswith('***For species') and not connMismatchFlag:
			connMismatchSpeciesRetry = connMismatchSpeciesRetry + 1
			connMismatchFlag=True
		    while '#' not in line:
			line = iin.readline()
			if line.startswith('***For species') and not connMismatchFlag:
			    connMismatchSpeciesRetry = connMismatchSpeciesRetry + 1
			    connMismatchFlag=True
			    #!note: I don't think this will catch a failure on the first attempt
		    attemptCounter=attemptCounter+1 #increment attempt# at each occurence of #
		if (line.startswith('Attempt')):
		    inchikey = line.split()[4]
		    if inchikey not in masterList: masterList.append(inchikey)
		    if attemptCounter > MopacKeywords :
			newGaussianArray[attemptCounter-MopacKeywords -2] = newGaussianArray[attemptCounter-MopacKeywords -2]+1 # note we subtract 2 from the index rather than 1 since when the program switches to Gaussian, it prints an extra spurious attempt #0
		    else:
			newMopacArray[attemptCounter-1] = newMopacArray[attemptCounter-1] + 1
		else:#starts with '*****Final attempt'
		    newFail = newFail+1
		    if (connMismatchFlag):
			connMismatchFail = connMismatchFail+1
		line = iin.readline()
	    elif (line.startswith('Pre-existing successful')):#a pre-existing successful read
		if (line.startswith('Pre-existing successful MOPAC')):
		    preexistMopac=preexistMopac+1
		    inchikey = line.split()[6]
		    if inchikey not in preexistMopacList: preexistMopacList.append(inchikey)
		    if inchikey not in masterList: masterList.append(inchikey)
		else:
		    preexistGaussian=preexistGaussian+1
		    inchikey = line.split()[5]
		    if inchikey not in preexistGaussianList: preexistGaussianList.append(inchikey)
		    if inchikey not in masterList: masterList.append(inchikey)
		line=iin.readline()
	    elif (line.startswith('Warning: Connectivity in quantum result for')):
		inchikey = line.split()[6]
		if inchikey not in connMismatchList: connMismatchList.append(inchikey)
		line=iin.readline()
	    else:
		line=iin.readline()
	iin.close()

	#compute totals


	#print the results
	print 'Number of species that completely failed QMTP = '+str(newFail)
	print 'Number of the complete failures that were due to apparent connectivity mismatches (CheckConnectivity: confirm) = '+ str(connMismatchFail)
	print 'Number of species that were retried due to apparent connectivity mismatches (CheckConnectivity: confirm) = '+str(connMismatchSpeciesRetry)
	print 'Number of reads of pre-existing successful MOPAC results = '+ str(preexistMopac)
	print 'Number of reads of pre-existing successful Gaussian results = '+ str(preexistGaussian)
	print 'Number of unique pre-existing MOPAC results read = '+ str(len(preexistMopacList))
	print 'Number of unique pre-existing Gaussian results read = '+str(len(preexistGaussianList))
	print 'Number of unique species with connectivity warnings (CheckConnectivity: check) = '+str(len(connMismatchList))
	print 'Number of species that finally succeeded at each attempt #:'
	newMopacTotal=0
	newGaussianTotal=0
	mopacTotalAttempts=newFail*MopacKeywords
	gaussianTotalAttempts=newFail*GaussianKeywords
	for i in range(len(newMopacArray)):
	    print 'Attempt #'+str(i+1)+' (MOPAC) : '+str(newMopacArray[i])
	    newMopacTotal = newMopacTotal + newMopacArray[i]
	    mopacTotalAttempts = mopacTotalAttempts + (i+1)*newMopacArray[i]
	for i in range(len(newGaussianArray)):
	    print 'Attempt #'+str(i+1)+' (Gaussian) : '+str(newGaussianArray[i])
	    newGaussianTotal = newGaussianTotal + newGaussianArray[i]
	    gaussianTotalAttempts = gaussianTotalAttempts + (i+1)*newGaussianArray[i]
	    mopacTotalAttempts = mopacTotalAttempts + MopacKeywords*newGaussianArray[i]
	print 'Number of species with new MOPAC results = ' + str(newMopacTotal)
	print 'Number of species with new Gaussian results = ' + str(newGaussianTotal)
	print 'Number of MOPAC attempts = ' + str(mopacTotalAttempts)
	print 'Number of Gaussian attempts = ' + str(gaussianTotalAttempts)
	print 'Total number of unique species with QMTP results (pre-existing and/or new) = ' + str(len(masterList))