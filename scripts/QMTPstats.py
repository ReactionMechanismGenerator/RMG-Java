# arg[1]: RMG log file from a run with QM with the "both" option
# arg[2]: 0 for "check"; 1 for "confirm"
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

#***For species CPGPQFYAGKHQTB-UHFFFAOYAZ an OpenBabel-based check suggests the optimized three-dimensional InChI (InChI=1/C5H6/c1-5-3-2-4-5/h2-4H2) does not match the intended (unmodified) InChI (InChI=1/C5H6/c1-2-4-5-3-1/h1-3H2). Will retry connectivity check with MoleCoor
#***For species CPGPQFYAGKHQTB-UHFFFAOYAZ a MoleCoor-based check suggests the optimized three-dimensional InChI (InChI=1/C5H6/c1-5-3-2-4-5/h2-4H2) does not match the intended (unmodified) InChI (InChI=1/C5H6/c1-2-4-5-3-1/h1-3H2).
#*****Attempt #1 on species CPGPQFYAGKHQTB-UHFFFAOYAZ (InChI=1/C5H6/c1-2-4-5-3-1/h1-3H2) failed. Will attempt a new keyword.
#***For species CPGPQFYAGKHQTB-UHFFFAOYAZ an OpenBabel-based check suggests the optimized three-dimensional InChI (InChI=1/C5H6/c1-5-3-2-4-5/h2-4H2) does not match the intended (unmodified) InChI (InChI=1/C5H6/c1-2-4-5-3-1/h1-3H2). Will retry connectivity check with MoleCoor
#***For species CPGPQFYAGKHQTB-UHFFFAOYAZ a MoleCoor-based check suggests the optimized three-dimensional InChI (InChI=1/C5H6/c1-5-3-2-4-5/h2-4H2) does not match the intended (unmodified) InChI (InChI=1/C5H6/c1-2-4-5-3-1/h1-3H2).
#*****Attempt #2 on species CPGPQFYAGKHQTB-UHFFFAOYAZ (InChI=1/C5H6/c1-2-4-5-3-1/h1-3H2) failed. Will attempt a new keyword.
#***For species CPGPQFYAGKHQTB-UHFFFAOYAZ an OpenBabel-based check suggests the optimized three-dimensional InChI (InChI=1/C5H6/c1-5-3-2-4-5/h2-4H2) does not match the intended (unmodified) InChI (InChI=1/C5H6/c1-2-4-5-3-1/h1-3H2). Will retry connectivity check with MoleCoor
#***For species CPGPQFYAGKHQTB-UHFFFAOYAZ a MoleCoor-based check suggests the optimized three-dimensional InChI (InChI=1/C5H6/c1-5-3-2-4-5/h2-4H2) does not match the intended (unmodified) InChI (InChI=1/C5H6/c1-2-4-5-3-1/h1-3H2).
#*****Attempt #3 on species CPGPQFYAGKHQTB-UHFFFAOYAZ (InChI=1/C5H6/c1-2-4-5-3-1/h1-3H2) failed. Will attempt a new keyword.
#Attempt #4 on species CPGPQFYAGKHQTB-UHFFFAOYAZ (InChI=1/C5H6/c1-2-4-5-3-1/h1-3H2) succeeded.
#Point group: C2v
#Thermo for CPGPQFYAGKHQTB-UHFFFAOYAZ: 127.6     67.2    18.37   24.32   29.54   33.86   40.39   45.06   52.1
	MopacKeywords = 10
	GaussianKeywords = 36
	if(int(sys.argv[2])==1):#0 for check, 1 for confirm
	    confirm = True
	else:
	    confirm = False

	#initialize counters
	preexistMopac=0
	preexistGaussian=0
	preexistMopacList=[]
	preexistGaussianList=[]
	masterList=[]
	connMismatchList=[]
	connMismatchBackupFixList=[]
	newMopacArray=[0 for x in range(MopacKeywords)]
	newGaussianArray=[0 for x in range(GaussianKeywords)]
	newMopacConnFixArray=[0 for x in range(MopacKeywords)]
	newGaussianConnFixArray=[0 for x in range(GaussianKeywords)]
	newFail=0
	connMismatchSpeciesRetry=0
	connMismatchFail=0
	connMismatchBackupFixesEarlierFail=0
	connMismatchBackupFixesOnFirstPrimaryFail=0

	#read in the file
	iin=open(sys.argv[1], 'r')
	line = iin.readline()
	while line != '':#'' marks end of file
	    #note: "***For species" should only be printed with check=confirm, so we don't need to check for that flag
	    if ('Attempt #1' in line or (line.startswith('***For species'))):#the beginning of a block for a new attempt
		#print line
		connMismatchFlag=False
		#dontreadmore=False
		#begin block (****used three places****)
		if (line.startswith('***For species')):
		    line = iin.readline()#read the next line to see whether the backup test failed also (in which case the next line would also start with '***For species'); otherwise, this should get to the attempt line or Warning line in the case of check
		    if (not line.startswith('***For species')):#the case was fixed by the backup method
			inchikey = line.split()[2]
			if inchikey not in connMismatchBackupFixList:
			    connMismatchBackupFixList.append(inchikey)
			    if connMismatchFlag:
				connMismatchBackupFixesEarlierFail=connMismatchBackupFixesEarlierFail+1
			    else:
				connMismatchBackupFixesOnFirstPrimaryFail=connMismatchBackupFixesOnFirstPrimaryFail+1
			line = iin.readline()#read the next line...should be 'For species' (without ***)
		    else:#else, the case was a genuine mismatch
			if (not connMismatchFlag):
			    connMismatchSpeciesRetry = connMismatchSpeciesRetry + 1
			    connMismatchFlag=True
			if not confirm:#read the next two lines to get to the Attempt line, (there will be a warning line)
			    line = iin.readline()
			    #at the this point, line should contain the warning
			    if not line.startswith('Warning: Connectivity in quantum result for'):
				print 'Algorithm error 0a:'+line
			    inchikey = line.split()[6]
			    if inchikey not in connMismatchList: connMismatchList.append(inchikey)
			    line = iin.readline()
			else:#read the next line to get to the Attempt line
			    line = iin.readline()
		#end block
		attemptCounter=1
		#at the end of here, line should contain "Attempt #1" (for confirm) (but only if there were no old (unconfirmed) files in QMfiles at the start of the run) or Pre-existing... or Attempt #1 (for check)
		if confirm:
		    #read until the line contains "pre-existing successful" or "Attempt #1"
		    while 'Attempt #1' not in line and not line.startswith('Pre-existing successful'):
			line = iin.readline()
		else:
		    if 'Attempt #1' not in line and not line.startswith('Pre-existing successful'):
			print 'Algorithm error 1:'+ line
		if (line.startswith('Pre-existing successful')):#a pre-existing successful read (following an (at least partial) mismatch
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
		else:#the Attempt #1 case
		    while (not line.startswith('Attempt') and not line.startswith('*****Final attempt') and not line == ''):#this loop gets to the successful attempt or the last attempt, counting the number of intermediate failed attempts in the process
			line = iin.readline()
			#begin block
			if (line.startswith('***For species')):
			    line = iin.readline()#read the next line to see whether the backup test failed also (in which case the next line would also start with '***For species'); otherwise, this should get to the attempt line or Warning line in the case of check
			    if (not line.startswith('***For species')):#the case was fixed by the backup method
				inchikey = line.split()[2]
				if inchikey not in connMismatchBackupFixList:
				    connMismatchBackupFixList.append(inchikey)
				    if connMismatchFlag:
					connMismatchBackupFixesEarlierFail=connMismatchBackupFixesEarlierFail+1
				    else:
					connMismatchBackupFixesOnFirstPrimaryFail=connMismatchBackupFixesOnFirstPrimaryFail+1
				line = iin.readline()#read the next line...should be 'For species' (without ***)
			    else:#else, the case was a genuine mismatch
				if (not connMismatchFlag):
				    connMismatchSpeciesRetry = connMismatchSpeciesRetry + 1
				    connMismatchFlag=True
				if not confirm:#read the next two lines to get to the Attempt line, (there will be a warning line)
				    line = iin.readline()
				    #at the this point, line should contain the warning
				    if not line.startswith('Warning: Connectivity in quantum result for'):
					print 'Algorithm error 0b:'+line
				    inchikey = line.split()[6]
				    if inchikey not in connMismatchList: connMismatchList.append(inchikey)
				    line = iin.readline()
				else:#read the next line to get to the Attempt line
				    line = iin.readline()
			    #at the this point, line should contain the next attempt"#"
			    if '#' not in line:
				print 'Algorithm error 2:'+line
			#end block
			while '#' not in line and line != '':#get to the next attempt line or "for species" line, ignoring, for example, imaginary frequencies warnings
			    line = iin.readline()
			    #begin block
			    if (line.startswith('***For species')):
				line = iin.readline()#read the next line to see whether the backup test failed also (in which case the next line would also start with '***For species'); otherwise, this should get to the attempt line or Warning line in the case of check
				if (not line.startswith('***For species')):#the case was fixed by the backup method
				    inchikey = line.split()[2]
				    if inchikey not in connMismatchBackupFixList:
					connMismatchBackupFixList.append(inchikey)
					if connMismatchFlag:
					    connMismatchBackupFixesEarlierFail=connMismatchBackupFixesEarlierFail+1
					else:
					    connMismatchBackupFixesOnFirstPrimaryFail=connMismatchBackupFixesOnFirstPrimaryFail+1
				    line = iin.readline()#read the next line...should be 'For species' (without ***)
				else:#else, the case was a genuine mismatch
				    if (not connMismatchFlag):
					connMismatchSpeciesRetry = connMismatchSpeciesRetry + 1
					connMismatchFlag=True
				    if not confirm:#read the next two lines to get to the Attempt line, (there will be a warning line)
					line = iin.readline()
					#at the this point, line should contain the warning
					if not line.startswith('Warning: Connectivity in quantum result for'):
					    print 'Algorithm error 0c:'+line
					inchikey = line.split()[6]
					if inchikey not in connMismatchList: connMismatchList.append(inchikey)
					line = iin.readline()
				    else:#read the next line to get to the Attempt line
					line = iin.readline()
				#at the this point, line should contain the next attempt"#"
				if '#' not in line:
				    print 'Algorithm error 3:'+line
			    #end block
			attemptCounter=attemptCounter+1 #increment attempt# at each occurence of #
		    if (line.startswith('Attempt')):
			inchikey = line.split()[4]
			if inchikey not in masterList: masterList.append(inchikey)
			if attemptCounter > MopacKeywords :
			    newGaussianArray[attemptCounter-MopacKeywords -2] = newGaussianArray[attemptCounter-MopacKeywords -2]+1 # note we subtract 2 from the index rather than 1 since when the program switches to Gaussian, it prints an extra spurious attempt #0
			    if connMismatchFlag:
				newGaussianConnFixArray[attemptCounter-MopacKeywords -2] = newGaussianConnFixArray[attemptCounter-MopacKeywords -2]+1
			else:
			    newMopacArray[attemptCounter-1] = newMopacArray[attemptCounter-1] + 1
			    if connMismatchFlag:
				newMopacConnFixArray[attemptCounter-1] = newMopacConnFixArray[attemptCounter-1]+1
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
	    else:
		line=iin.readline()
	iin.close()

	#print the results
	print 'Number of species that completely failed QMTP = '+str(newFail)
	if confirm: print 'Number of the complete failures that were due to apparent connectivity mismatches (CheckConnectivity: confirm) = '+ str(connMismatchFail)
	if confirm: print 'Number of species that were retried due to apparent connectivity mismatches (CheckConnectivity: confirm) = '+str(connMismatchSpeciesRetry)
	if not confirm: print 'Number of unique species with connectivity warnings (CheckConnectivity: check) = '+str(len(connMismatchList))
	print 'Number of cases where backup connectivity check identified a match missed by the primary method (early/late) = '+str(len(connMismatchBackupFixList))+' ('+str(connMismatchBackupFixesOnFirstPrimaryFail)+'/'+str(connMismatchBackupFixesEarlierFail)+')'#this could double-count for "check" case
	print 'Number of reads of pre-existing successful MOPAC results = '+ str(preexistMopac)
	print 'Number of reads of pre-existing successful Gaussian results = '+ str(preexistGaussian)
	print 'Number of unique pre-existing MOPAC results read = '+ str(len(preexistMopacList))
	print 'Number of unique pre-existing Gaussian results read = '+str(len(preexistGaussianList))
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
	if confirm:
	    print 'Number of species with connectivity apparently fixed at each attempt #:'
	    for i in range(len(newMopacConnFixArray)):
		print 'Attempt #'+str(i+1)+' (MOPAC) : '+str(newMopacConnFixArray[i])
	    for i in range(len(newGaussianConnFixArray)):
		print 'Attempt #'+str(i+1)+' (Gaussian) : '+str(newGaussianConnFixArray[i])
	print 'Number of species with new MOPAC results = ' + str(newMopacTotal)
	print 'Number of species with new Gaussian results = ' + str(newGaussianTotal)
	print 'Number of MOPAC attempts = ' + str(mopacTotalAttempts)
	print 'Number of Gaussian attempts = ' + str(gaussianTotalAttempts)
	print 'Total number of unique species with QMTP results (pre-existing and/or new) = ' + str(len(masterList))