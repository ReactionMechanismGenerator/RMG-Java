/*
 * FrequencyGroups.java 
 * created by Greg Magoon 11/17/08 using GATP.java as a starting point
 * This class is the frequency database analogue to GATP
 */

package jing.chem;



import java.util.*;
import jing.chemUtil.*;


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
   //FreqData class temporarily commented out for testing 11/18/08
    // public FreqData generateFreqData(ChemGraph p_chemGraph) {
      public LinkedList generateFreqData(ChemGraph p_chemGraph){
       // result = new FreqData();//note: FreqData class has not yet been created, but may be something similar to ThreeFrequencyModel
        LinkedList groupCount=getFreqGroup(p_chemGraph);//use database to count groups in structure
        //write input to Franklin's code (will also need to have heat capacity, internal rotor info, etc. passed to this function) ****not yet done
        //call Franklin's code****not yet done
        //read in results of Franklin's code****not yet done
       // return result;
          return groupCount;
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

        //from groupCountMap, create the "result" LinkedList consisting of the numbers of each type of group used by Franklin's code in the order that his code requires them
        //if the group name is not in groupCountMap, a zero will be used
        //in the future, we may want to "un-hardcode" this by also reading in a file with a list of the different groups used by Franklin's code and the order in which they occur
        //*this section of code not yet done*
        
         p_chemGraph.setCentralNode(oldCentralNode);
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
