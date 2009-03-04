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
            //// the following line won't compile, with error:
            ////  generateFreqData(jing.chem.ChemGraph,jing.chem.ThermoData) in jing.chem.FrequencyGroups cannot be applied to (jing.chem.ChemGraph)
            //// rwest commenting out from the CVS version as it's blocking my automated build process
            // freqGroups.generateFreqData(cg2);	
	  			
	    	}
	    	catch (IOException e){
	    		System.out.println(e.toString());
	    	}

	    	
	    }
}
