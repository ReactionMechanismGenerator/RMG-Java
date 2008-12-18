/*
 * FrequencyGroups.java 
 * created by Greg Magoon 11/17/08 using GATP.java as a starting point
 * This class is the frequency database analogue to GATP
 */

package jing.chem;



import java.util.*;
import jing.chemUtil.*;

import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.BufferedReader;

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
		int degeneracy = 0;
		
		degeneracy = 
				9*((Integer) groupCount.get(0)) + 5*((Integer) groupCount.get(1)) + 
				3*((Integer) groupCount.get(2)) + 7*((Integer) groupCount.get(3)) + 
				5*((Integer) groupCount.get(4))	+ 5*((Integer) groupCount.get(5)) + 
				3*((Integer) groupCount.get(6)) + 3*((Integer) groupCount.get(7)) + 
				2*((Integer) groupCount.get(8)) + 6*((Integer) groupCount.get(9)) + 
				4*((Integer) groupCount.get(10)) + 4*((Integer) groupCount.get(11)) + 
				5*((Integer) groupCount.get(12)) + 6*((Integer) groupCount.get(13)) + 
                3*((Integer) groupCount.get(14)) + 4*((Integer) groupCount.get(15)) + 
				2*((Integer) groupCount.get(16)) + 2*((Integer) groupCount.get(17)) + 
				3*((Integer) groupCount.get(18)) + 4*((Integer) groupCount.get(19)) + 
				1*((Integer) groupCount.get(20)) + 5*((Integer) groupCount.get(21)) + 
				3*((Integer) groupCount.get(22)) + 3*((Integer) groupCount.get(23));
		
		
		int nFreq = 3 * atoms - 5 - rotor - degeneracy - linearity;
		if (nFreq < 0) {
			//System.out.println("ERROR: Attempted to fit a negative number of frequencies!");
			//System.exit(0);
			rotor = 0;
		}
		
		//(file writing code based on code in JDASSL.java)
        File franklInput = new File("frankie/dat");
        try {
            FileWriter fw = new FileWriter(franklInput);
            fw.write(p_thermoData.getCp300()+" "+p_thermoData.getCp400()+" "+p_thermoData.getCp500()+" "+p_thermoData.getCp600()+" "+p_thermoData.getCp800()+" "+p_thermoData.getCp1000()+" "+p_thermoData.getCp1500()+"\n");
            fw.write(atoms+"\n");
            fw.write(rotor+"\n");
            fw.write(linearity+"\n");
            //print the group counts to the file
            for (Iterator iter = groupCount.iterator(); iter.hasNext(); ) {
				Integer i = (Integer) iter.next();
				fw.write(i + "\n");
			}
			fw.close();
        }
        catch (IOException e) {
            System.err.println("Problem writing frequency estimation input file!");
            e.printStackTrace();
        }
		
		touchOutputFile();
		
        //call Franklin's code
        try{
            String dir = System.getProperty("RMG.workingDirectory");
			File runningdir=new File("frankie/");
            String command = dir + "/software/frankie/frankie.exe";
           	Process freqProc = Runtime.getRuntime().exec(command, null, runningdir); 
            int exitValue = freqProc.waitFor();
        }
        catch (Exception e) {
            String err = "Error in running frequency estimation process \n";
            err += e.toString();
            e.printStackTrace();
        }
        
         //read in results of Franklin's code into result
        File franklOutput = new File("frankie/rho_input");
        String line="";
        try { 
            
			FileReader fr = new FileReader(franklOutput);
            BufferedReader br = new BufferedReader(fr);
            
			// Read information about molecule (numuber of atoms, number of rotors, linearity)
			atoms = Integer.parseInt(br.readLine().trim());
            int nHind = Integer.parseInt(br.readLine().trim());
            int nonlinearityQ = Integer.parseInt(br.readLine().trim());
            
			// Determine the expected number of harmonic oscillator frequencies
			nFreq = 0;
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
					i++;
				}
			}
			if (i != nFreq) {
				System.out.println("Warning: Number of frequencies read is less than expected.");
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
					i++;
				}
			}
			if (i != nHind) {
				System.out.println("Warning: Number of hindered rotors read is less than expected.");
			}
			
			fr.close();
        }
        catch (IOException e) {
                System.err.println("Problem reading frequency estimation output file!");
                e.printStackTrace();           
        }
		catch (NullPointerException e) {
                System.err.println("Problem reading frequency estimation output file!");
                e.printStackTrace();           
        }
        
		// Rename input and output files
		String newName = species.getName()+"("+String.valueOf(species.getID())+")";
        File f = new File("frankie/dat");
		File newFile = new File("frankie/"+newName+"_input");
		f.renameTo(newFile);
		f = new File("frankie/rho_input");
		newFile = new File("frankie/"+newName+"_output");
		f.renameTo(newFile);
		
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


    protected static FrequencyGroups getINSTANCE() {
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
			System.exit(0);
		}
		catch(Exception e) {
			System.out.println(e.getMessage());
		}
	}


}
/*********************************************************************
	File Path	: RMG\RMG\jing\chem\FrequencyGroups.java
*********************************************************************/
