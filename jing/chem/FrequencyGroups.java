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
      public SpectroscopicData generateFreqData(ChemGraph p_chemGraph, ThermoData p_thermoData) {
        LinkedList groupCount=getFreqGroup(p_chemGraph);//use database to count groups in structure
        double[] rotfreq=null;
        double[] rotV=null;
        double[] vibs=null;
        
        //write input to Franklin's code
        //determine linearity by mapping true/false to 0/1
        int linearity=1;
        if(p_chemGraph.isLinear())
            linearity=0;
        //(file writing code based on code in JDASSL.java)
        File franklInput = new File("dat");
        try {
            FileWriter fw = new FileWriter(franklInput);
            fw.write(p_thermoData.getCp300()+" "+p_thermoData.getCp400()+" "+p_thermoData.getCp500()+" "+p_thermoData.getCp600()+" "+p_thermoData.getCp800()+" "+p_thermoData.getCp1000()+" "+p_thermoData.getCp1500()+"\n");
            fw.write(p_chemGraph.getAtomNumber()+"\n");
            fw.write(p_chemGraph.getInternalRotor()+"\n");
            fw.write(linearity+"\n");
            //print the group counts to the file
            Iterator iter = groupCount.iterator();
            while (iter.hasNext())
                fw.write((Integer)iter.next());
            fw.close();
        }
        catch (IOException e) {
            System.err.println("Problem writing frequency estimation input file!");
            e.printStackTrace();
        }
		
        //call Franklin's code
        try{
            File runningdir=new File("?????");
            String command = "?????";
            Process freqProc = Runtime.getRuntime().exec(command, null, runningdir); 
            int exitValue = freqProc.waitFor();
        }
        catch (Exception e) {
            String err = "Error in running frequency estimation process \n";
            err += e.toString();
            e.printStackTrace();
        }
        
         //read in results of Franklin's code into result
        File franklOutput = new File("rho_input");
        String line="";
        try{
            FileReader fr = new FileReader(franklOutput);
            BufferedReader br = new BufferedReader(fr);
            line = br.readLine();
            int atoms=Integer.parseInt(line.trim());
            line = br.readLine();
            int rotors=Integer.parseInt(line.trim());
            line = br.readLine();
            int nonlinearityQ=Integer.parseInt(line.trim());
            //read the HO frequencies, RR frequencies, and barrier heights
            int nfreq=0;
            if (nonlinearityQ==1)
                nfreq=3*atoms-6;
            else
                nfreq=3*atoms-5;
            vibs=new double[nfreq];
            line=br.readLine();
            int i=0;
            StringTokenizer st=new StringTokenizer(line.trim());
            while(st.hasMoreTokens()){
                vibs[i]=Double.parseDouble(st.nextToken());
                i++;
            }
            //at this point, i should equal nfreqs; we could put a check here to confirm this
            rotfreq=new double[rotors];
            rotV=new double[rotors];
            line=br.readLine();
            i=0;
            StringTokenizer st2=new StringTokenizer(line.trim());
            while(st2.hasMoreTokens()){
                rotfreq[i]=Double.parseDouble(st2.nextToken());
                rotV[i]=Double.parseDouble(st2.nextToken());
                i++;
            }
            //at this point, i should equal rotors; we could put a check here to confirm this
            fr.close();
        }
        catch (IOException e) {
                System.err.println("Problem reading frequency estimation output file!");
                e.printStackTrace();           
        }
        
        SpectroscopicData result = new SpectroscopicData(vibs, rotfreq, rotV);
        
        return result;
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


}
/*********************************************************************
	File Path	: RMG\RMG\jing\chem\FrequencyGroups.java
*********************************************************************/
