/*
 * gmagoon 11/18/08: tests FrequencyGroups.java...based on Sandeep's IdentifyMatchesSites.java
 */

package jing.chem;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import jing.chemParser.ChemParser;
import jing.chemParser.InvalidGraphFormatException;
import jing.chemUtil.Graph;


public class FrequencyGroupTester {
	 public static void main(String[] args) {
	    	File funcGroups = new File ("c:/Users/User1/Documents/NetBeansProjects/RMGdev/build/databases/RMG_database/frequencies/testmolecule.txt");
	    	try{
	    		BufferedReader reader = new BufferedReader(new FileReader(funcGroups));
	    		Graph Graph2 = null;
                        try {
                                Graph2 = ChemParser.readChemGraph(reader);
                        }
                        catch (InvalidGraphFormatException e) {
                                throw new InvalidFunctionalGroupException(e.getMessage());
                        }
                        ChemGraph cg2 = null;
                        try {
                                cg2 = ChemGraph.make(Graph2);
                        } catch (InvalidChemGraphException e) {
                                // TODO Auto-generated catch block
                                e.printStackTrace();
                        } catch (ForbiddenStructureException e) {
                                // TODO Auto-generated catch block
                                e.printStackTrace();
                        }
	  		FrequencyGroups freqGroups= FrequencyGroups.getINSTANCE();
            //// the following line won't compile
            ////  generateFreqData(jing.chem.ChemGraph,jing.chem.ThermoData) in jing.chem.FrequencyGroups cannot be applied to (jing.chem.ChemGraph)
            //// rwest commenting out from the CVS version as it's blocking my automated build process
            // freqGroups.generateFreqData(cg2);	
	  			
	    	}
	    	catch (IOException e){
	    		System.out.println(e.toString());
	    	}

	    	
	    }
}
