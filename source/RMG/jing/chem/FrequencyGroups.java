////////////////////////////////////////////////////////////////////////////////
//
//	RMG - Reaction Mechanism Generator
//
//	Copyright (c) 2002-2009 Prof. William H. Green (whgreen@mit.edu) and the
//	RMG Team (rmg_dev@mit.edu)
//
//	Permission is hereby granted, free of charge, to any person obtaining a
//	copy of this software and associated documentation files (the "Software"),
//	to deal in the Software without restriction, including without limitation
//	the rights to use, copy, modify, merge, publish, distribute, sublicense,
//	and/or sell copies of the Software, and to permit persons to whom the
//	Software is furnished to do so, subject to the following conditions:
//
//	The above copyright notice and this permission notice shall be included in
//	all copies or substantial portions of the Software.
//
//	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//	FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//	DEALINGS IN THE SOFTWARE.
//
////////////////////////////////////////////////////////////////////////////////



package jing.chem;



import java.util.*;
import jing.chemUtil.*;

import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;

public class FrequencyGroups{//gmagoon 111708: removed "implements GeneralGAPP"

    private static FrequencyGroups INSTANCE = new FrequencyGroups();		//## attribute INSTANCE
    protected FrequencyDatabase freqLibrary;

    // Constructors

    //## operation FrequencyGroups()
    private FrequencyGroups() {
        initFreqLibrary();
        //#]
    }

   //gmagoon 11/17/08: based off of generateThermoData from GATP.java;
   // this function generates a list of frequencies using Franklin's code for use in pressure-dependent network calculations
      public SpectroscopicData generateFreqData(Species species) {
		
		ChemGraph p_chemGraph = species.getChemGraph();
		ThermoData p_thermoData = species.getThermoData();
				
		// Skip if species is monatomic
		if (p_chemGraph.getAtomNumber() < 2)
			return new SpectroscopicData();
		  
		LinkedList groupCount=getFreqGroup(p_chemGraph);//use database to count groups in structure
        
		double[] vibFreq = null;
        double[] hindFreq = null;
		double[] hindBarrier = null;
			
		int atoms = p_chemGraph.getAtomNumber();
		int rotor = p_chemGraph.getInternalRotor();
		int linearity = (p_chemGraph.isLinear()) ? 0 : 1;	// 0 if linear, 1 if nonlinear
                if(p_chemGraph.isAcyclic()){//gmagoon 2/20/09: only perform check if molecule is acyclic; for cyclic cases, groups are still counted and constraints on degeneracy will not necessarily be enforced
                    int degeneracy = getDegeneracy(groupCount);

                    int nFreq = 3 * atoms - 5 - rotor - degeneracy - linearity;
                    if (nFreq < 0) {
                            System.out.println(species.getName() + "(" + Integer.toString(species.getID()) +
                                    ") is overspecified: " +
                                            Integer.toString(degeneracy) + " harmonic oscillators and " +
                                            Integer.toString(rotor) + " internal rotors are specified, but only " +
                                            Integer.toString(3 * atoms - 5 - linearity) + " modes are allowed.");
                            // If possible, turn off internal rotors until the number of fitted frequencies is zero
                            if (nFreq + rotor >= 0) {
                                    System.out.println("Turning off " + Integer.toString(Math.abs(nFreq)) + 
                                                    " internal rotors.");
                                    rotor += nFreq;
                            }
                            else {
                                    // Turn off functional groups until the system is underspecified
                                    System.out.println("Turning off functional groups to make problem underspecified.");
                                    while (nFreq < 0) {
                                            removeFunctionalGroup(groupCount);
                                            degeneracy = getDegeneracy(groupCount);
                                            nFreq = 3 * atoms - 5 - rotor - degeneracy - linearity;
                                    }
                            }
                    }
                 }
		
		//(file writing code based on code in JDASSL.java)
        File franklInput = new File("frankie/dat");
        try {
            FileWriter fw = new FileWriter(franklInput);
            fw.write(p_thermoData.getCp300()+" "+p_thermoData.getCp400()+" "+p_thermoData.getCp500()+" "+p_thermoData.getCp600()+" "+p_thermoData.getCp800()+" "+p_thermoData.getCp1000()+" "+p_thermoData.getCp1500()+"\n");
            fw.write(atoms+"\n");
//            fw.write(rotor+"\n");
            // Commented out by MRH (in conjunction with CFG) on 12-Jun-2009
            //	Some cyclic species were taking ~30 seconds to converge
            //	CFG suggested setting the rotor number to zero for all cyclic species
            if (p_chemGraph.isAcyclic()) fw.write(rotor+"\n");
            else fw.write(0+"\n");
            fw.write(linearity+"\n");
            //print the group counts to the file for acyclic case
            if(p_chemGraph.isAcyclic()){
                for (Iterator iter = groupCount.iterator(); iter.hasNext(); ) {
                    Integer i = (Integer) iter.next();
                    fw.write(i + "\n");
		}
                fw.write(0+"\n");//write a zero for acyclic species
            }
            else{//cyclic case: print zeroes for group values instead, followed by a line with the number of hydrogens
                for (Iterator iter = groupCount.iterator(); iter.hasNext(); ) {
                    Integer i = (Integer) iter.next();
                    fw.write(0 + "\n");
                }
                fw.write(p_chemGraph.getHydrogenNumber()+"\n");
            }
           
                
            fw.close();

        }
        catch (IOException e) {
            System.err.println("Problem writing frequency estimation input file!");
            e.printStackTrace();
        }
		
		touchOutputFile();
		  boolean frankieSuccess = false;
		  int frankieOutputFlag = 0;
		//call Franklin's code
        try{
            String dir = System.getProperty("RMG.workingDirectory");
            File runningdir=new File("frankie");
            String[] command = { dir + "/bin/frankie.exe" };
            Process freqProc = Runtime.getRuntime().exec(command, null, runningdir); 
            InputStream is = freqProc.getInputStream();
            InputStreamReader isr = new InputStreamReader(is);
            BufferedReader br = new BufferedReader(isr);
            String line=null;
            while ( (line = br.readLine()) != null) {
				line = line.trim();
				if (line.contains("Output flag from DQED, IGO ="))
					frankieOutputFlag =  Integer.parseInt( line.substring(line.length()-1) );
            }
            int exitVal = freqProc.waitFor();
		
        }
        catch (Exception e) {
            String err = "Error in running frequency estimation process \n";
            err += e.toString();
            e.printStackTrace();
        }
		  
		  if (frankieOutputFlag == 4) 
			  frankieSuccess = true;
		  if (frankieOutputFlag == 8) 
			  System.err.println("Frankie exceeded maximum number of iterations");	  
		  
		  if (!frankieSuccess) {
			  System.err.println("Frankie.exe wasn't fully successful: "+ String.format("species %2$d had output flag %1$d",frankieOutputFlag,species.getID() ));
			  System.err.println(String.format("Saving input file as 'frankie/dat.%d.%d' should you wish to debug.",frankieOutputFlag,species.getID() ));	
			  franklInput.renameTo( new File(String.format("frankie/dat.%d.%d",frankieOutputFlag,species.getID() )) );			  
		  }

        
         //read in results of Franklin's code into result
        File franklOutput = new File("frankie/rho_input");
        String line="";
        try { 
            
			FileReader fr = new FileReader(franklOutput);
            BufferedReader br = new BufferedReader(fr);
			
			// If the output file is empty, then an error occurred
			if (franklOutput.length() == 0)
				throw new IOException("Error: FrequencyGroups estimation for " + 
						species.getName() + "(" + Integer.toString(species.getID()) + ") failed.");
            
			// Read information about molecule (numuber of atoms, number of rotors, linearity)
			atoms = Integer.parseInt(br.readLine().trim());
            int nHind = Integer.parseInt(br.readLine().trim());
            int nonlinearityQ = Integer.parseInt(br.readLine().trim());
            
			// Determine the expected number of harmonic oscillator frequencies
			int nFreq = 0;
            if (nonlinearityQ == 1)
                nFreq = 3 * atoms - 6 - nHind;
            else
                nFreq = 3 * atoms - 5 - nHind;
			
			// Read the harmonic oscillator frequencies
			vibFreq = new double[nFreq];
            int i = 0;
			while (i < nFreq) {
				StringTokenizer st = new StringTokenizer(br.readLine().trim());
				if (st.countTokens() == 0)
					continue;
				else {
					vibFreq[i] = Double.parseDouble(st.nextToken());
					if (vibFreq[i] < 0)
						throw new IOException("Encountered a negative vibrational frequency while reading Frankie output.");
					i++;
				}
			}
			if (i != nFreq) {
				throw new IOException("Number of frequencies read is less than expected.");
			}
			
			// Read the hindered rotor frequencies and barrier heights
			hindFreq = new double[nHind];
			hindBarrier = new double[nHind];
			i = 0;
			while (i < nHind) {
				StringTokenizer st = new StringTokenizer(br.readLine().trim());
				if (st.countTokens() == 0)
					continue;
				else {
					hindFreq[i] = Double.parseDouble(st.nextToken());
					hindBarrier[i] = Double.parseDouble(st.nextToken());
					if (hindFreq[i] < 0)
						throw new IOException("Encountered a negative hindered rotor frequency while reading Frankie output.");
					if (hindBarrier[i] < 0)
						throw new IOException("Encountered a negative hindered rotor barrier while reading Frankie output.");
					i++;
				}
			}
			if (i != nHind) {
				throw new IOException("Number of hindered rotors read is less than expected.");
			}
			
			fr.close();
        }
        catch (IOException e) {
                System.err.println("Problem reading frequency estimation output file!");
				System.out.println(e.getMessage());
                e.printStackTrace();
				System.exit(0);
        }
		catch (NullPointerException e) {
                System.err.println("Problem reading frequency estimation output file!");
                e.printStackTrace();           
        }
        
		// Rename input and output files
		/*String newName = species.getName()+"("+String.valueOf(species.getID())+")";
        File f = new File("frankie/dat");
		File newFile = new File("frankie/"+newName+"_input");
		f.renameTo(newFile);
		f = new File("frankie/rho_input");
		newFile = new File("frankie/"+newName+"_output");
		f.renameTo(newFile);*/
		
        return new SpectroscopicData(vibFreq, hindFreq, hindBarrier);
        //#]
    }

    //11/17/08 gmagoon: modified from getGAGroup in class GATP
    //this function will likely need to be modified frequently to adapt to changes in Franklin's code
    /**
    Requires: pass-in ChemGraph object repOk() == true;
    Effects: counts the number of groups in the order required by Franklin's code
    Modifies:
    */
    public LinkedList getFreqGroup(ChemGraph p_chemGraph) {
        LinkedList result = new LinkedList();
        HashMap oldCentralNode = (HashMap)(p_chemGraph.getCentralNode()).clone();

        // find all the groups, and store them in the HashMap groupCountMap, where the name of the group is the key and the value is the number of times the group has been encountered in the structure
        HashMap groupCountMap = new HashMap();
        Iterator iter = p_chemGraph.getNodeList();
        while (iter.hasNext()) {
          	Node node = (Node)iter.next();
          	Atom atom = (Atom)node.getElement();
          	if (!(atom.getType().equals("H"))) {
                    p_chemGraph.resetThermoSite(node);
                    String thisFreqGroupName = freqLibrary.findFreqGroupName(p_chemGraph);
                    //if groupCountMap already contains the group name, increment the value by one; otherwise, add the group name as a key with an initial value of one
                    if(groupCountMap.containsKey(thisFreqGroupName)){
                        Integer oldCount=(Integer)(groupCountMap.get(thisFreqGroupName));
                        groupCountMap.put(thisFreqGroupName, oldCount+1);
                    }
                    else{
                        groupCountMap.put(thisFreqGroupName,1);
                    }
           	}
        }
        p_chemGraph.setCentralNode(oldCentralNode);
        
        //from groupCountMap, create the "result" LinkedList consisting of the numbers of each type of group used by Franklin's code in the order that his code requires them
        //if the group name is not in groupCountMap, a zero will be used
        //in the future, we may want to "un-hardcode" this by also reading in a file with a list of the different groups used by Franklin's code and the order in which they occur
        String[] orderedInputGroups={"RsCH3","RdCH2","CtCH","RsCH2sR","CdCHsR","Aldehyde","Cumulene","Ketene","CtCsR","RsCHsR2","CdCsR2","Ketone","RsCsR3","RsCH2r","RdCHr","RsCHrsR","CdCrsR","OdCrsR","RsCrsR2","Alcohol","Ether","ROOH","ROOR","Peroxy"};//this should contain the group names (or keys names) used by Franklin's frequency estimation code in the order that his input format requires them
        for(int i=1;i<=orderedInputGroups.length;i++){
            String inputGroup=orderedInputGroups[i-1];
            if(groupCountMap.containsKey(inputGroup))
                result.add((Integer)(groupCountMap.get(inputGroup)));
            else
                result.add(0);
        }
        
        return result;

        //#]
    }


    //11/17/08 gmagoon: modifed from initGAGroupLibrary from GATP.java
    protected void initFreqLibrary() {
        freqLibrary = FrequencyDatabase.getINSTANCE();
        //#]
    }


    public static FrequencyGroups getINSTANCE() {
        return INSTANCE;
    }
	
	private void touchOutputFile() {
		try {
			File output = new File("frankie/rho_input");
			if (output.exists())
				output.delete();
			output.createNewFile();
		}
		catch(IOException e) {
			System.out.println("Error: Unable to touch file \"frankie/rho_input\".");
			System.out.println(e.getMessage());
			System.exit(0);
		}
		catch(Exception e) {
			System.out.println(e.getMessage());
		}
	}

	private int getDegeneracy(LinkedList groupCount) {
		return	8*((Integer) groupCount.get(0)) + 5*((Integer) groupCount.get(1)) +
				3*((Integer) groupCount.get(2)) + 7*((Integer) groupCount.get(3)) +
				5*((Integer) groupCount.get(4))	+ 5*((Integer) groupCount.get(5)) +
				3*((Integer) groupCount.get(6)) + 3*((Integer) groupCount.get(7)) +
				2*((Integer) groupCount.get(8)) + 6*((Integer) groupCount.get(9)) +
				4*((Integer) groupCount.get(10)) + 4*((Integer) groupCount.get(11)) +
				5*((Integer) groupCount.get(12)) + 5*((Integer) groupCount.get(13)) +
                3*((Integer) groupCount.get(14)) + 4*((Integer) groupCount.get(15)) +
				2*((Integer) groupCount.get(16)) + 2*((Integer) groupCount.get(17)) +
				3*((Integer) groupCount.get(18)) + 2*((Integer) groupCount.get(19)) +
				1*((Integer) groupCount.get(20)) + 4*((Integer) groupCount.get(21)) +
				1*((Integer) groupCount.get(22)) + 2*((Integer) groupCount.get(23));//gmagoon: updated for Franklin's new code 2/20/09 (Hydrogen group for rings not included)
	}

	private void removeFunctionalGroup(LinkedList groupCount) {
		if ((Integer) groupCount.get(20) > 0) groupCount.set(20, (Integer) groupCount.get(20) - 1);
		else if ((Integer) groupCount.get(8) > 0) groupCount.set(8, (Integer) groupCount.get(8) - 1);
		else if ((Integer) groupCount.get(16) > 0) groupCount.set(16, (Integer) groupCount.get(16) - 1);
		else if ((Integer) groupCount.get(17) > 0) groupCount.set(17, (Integer) groupCount.get(17) - 1);
		else if ((Integer) groupCount.get(2) > 0) groupCount.set(2, (Integer) groupCount.get(2) - 1);
		else if ((Integer) groupCount.get(6) > 0) groupCount.set(6, (Integer) groupCount.get(6) - 1);
		else if ((Integer) groupCount.get(7) > 0) groupCount.set(7, (Integer) groupCount.get(7) - 1);
		else if ((Integer) groupCount.get(14) > 0) groupCount.set(14, (Integer) groupCount.get(14) - 1);
		else if ((Integer) groupCount.get(18) > 0) groupCount.set(18, (Integer) groupCount.get(18) - 1);
		else if ((Integer) groupCount.get(23) > 0) groupCount.set(23, (Integer) groupCount.get(23) - 1);
		else if ((Integer) groupCount.get(10) > 0) groupCount.set(10, (Integer) groupCount.get(10) - 1);
		else if ((Integer) groupCount.get(11) > 0) groupCount.set(11, (Integer) groupCount.get(11) - 1);
		else if ((Integer) groupCount.get(15) > 0) groupCount.set(15, (Integer) groupCount.get(15) - 1);
		else if ((Integer) groupCount.get(19) > 0) groupCount.set(19, (Integer) groupCount.get(19) - 1);
		else if ((Integer) groupCount.get(1) > 0) groupCount.set(1, (Integer) groupCount.get(1) - 1);
		else if ((Integer) groupCount.get(4) > 0) groupCount.set(4, (Integer) groupCount.get(4) - 1);
		else if ((Integer) groupCount.get(5) > 0) groupCount.set(5, (Integer) groupCount.get(5) - 1);
		else if ((Integer) groupCount.get(12) > 0) groupCount.set(12, (Integer) groupCount.get(12) - 1);
		else if ((Integer) groupCount.get(21) > 0) groupCount.set(21, (Integer) groupCount.get(21) - 1);
		else if ((Integer) groupCount.get(9) > 0) groupCount.set(9, (Integer) groupCount.get(9) - 1);
		else if ((Integer) groupCount.get(13) > 0) groupCount.set(13, (Integer) groupCount.get(13) - 1);
		else if ((Integer) groupCount.get(3) > 0) groupCount.set(3, (Integer) groupCount.get(3) - 1);
		else if ((Integer) groupCount.get(0) > 0) groupCount.set(0, (Integer) groupCount.get(0) - 1);
	}

}
/*********************************************************************
	File Path	: RMG\RMG\jing\chem\FrequencyGroups.java
*********************************************************************/
