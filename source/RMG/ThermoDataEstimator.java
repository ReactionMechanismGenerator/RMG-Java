////////////////////////////////////////////////////////////////////////////////
//
//	RMG - Reaction Mechanism Generator
//
//	Copyright (c) 2002-2011 Prof. William H. Green (whgreen@mit.edu) and the
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

import java.util.*;
import java.io.*;

import jing.chem.*;
import jing.chemParser.*;
import jing.param.*;
import jing.chemUtil.*;
//import bondGroups.*;
import jing.rxnSys.*;


public class ThermoDataEstimator {//gmagoon 7/24/09: based off of Thermo.java rev. 1.6; this performs original functionality described in the manual

  //## configuration RMG::RMG
  //first argument will be file to read
public static void main(String[] args) {
//  initializeSystemProperties();
//
//          initializeSystemProperties();
	RMG.globalInitializeSystemProperties();
	
	File GATPFit = new File("GATPFit");
	GATPFit.mkdir();

		LinkedList<ChemGraph> graphList = new LinkedList<ChemGraph>();
                LinkedList<String> nameList = new LinkedList<String>();
                LinkedHashMap speciesFromInputFile = new LinkedHashMap();

		File file = new File(args[0]);

		try {			
			BufferedReader reader = new BufferedReader(new FileReader(file));

                        String line = ChemParser.readMeaningfulLine(reader, true);
                        if (line.toLowerCase().startsWith("database")) {
                            RMG.extractAndSetDatabasePath(line);
                        }
                        else {
                            System.err.println("ThermoDataEstimator: Could not"
                                    + " locate the Database field");
                            System.exit(0);
                        }

			/*
			 * 14MAR2010: Allowing ThermoDataEstimator to read from Primary Reaction Library
			 */
			ReactionModelGenerator rmg = new ReactionModelGenerator();
            line = ChemParser.readMeaningfulLine(reader, true);
            if (line.toLowerCase().startsWith("primarythermolibrary")) {
            	rmg.readAndMakePTL(reader);
            }
            else {
            	System.err.println("ThermoDataEstimator: Could not locate the PrimaryThermoLibrary field." +
            			"Line read was: " + line);
            	System.exit(0);
            }
			
			// Read adjacency lists from file until an exception is thrown
            line = ChemParser.readMeaningfulLine(reader, true);
			while (line != null) {
                            nameList.add(line);
                            Graph g = ChemParser.readChemGraph(reader);
				ChemGraph cg = ChemGraph.make(g);

                                ReactionModelGenerator.addChemGraphToListIfNotPresent_ElseTerminate(speciesFromInputFile,cg,"");

				graphList.add(cg);
				line = ChemParser.readMeaningfulLine(reader, true);
			}

		}
		catch (InvalidChemGraphException e) {
			e.printStackTrace();
		} catch (ForbiddenStructureException e) {
			e.printStackTrace();
		}
		catch (IOException e) {
			System.out.println(e.toString());
		}

                int counter = 0;
		for (ListIterator<ChemGraph> iter = graphList.listIterator(); iter.hasNext(); ) {
			ChemGraph chemgraph = iter.next();


         Species spe = Species.make(nameList.get(counter),chemgraph);
		 /*
		  *	Following line added by MRH on 10Aug2009:
		  *		After the species is made, the chemgraph is not necessarily the same
		  *		(e.g. the species contains resonance isomers, and the adjacency list
		  *		supplied by the user is not the most stable isomer).  This was causing
		  *		a discrepancy to be displayed to the screen: the "H=" value we were
		  *		displaying corresponded to the chemgraph supplied by the user; the 
		  *		ThermData corresponded to the most stable isomer.
		  *
		  *		An example of a troublesome chemgraph:
		  *		1 C 0 {2,T}
		  *		2 C 0 {1,T} {3,S}
		  *		3 O 0 {2,S}
		  */
		 chemgraph = spe.getChemGraph();
		 System.out.println(chemgraph);
         System.out.println("The number of resonance isomers is " + spe.getResonanceIsomersHashSet().size());
		 System.out.println("The NASA data is \n!"+spe.getNasaThermoSource()+"\n"+
                         "!" + chemgraph.getThermoComments() + "\n" +
				 spe.getNasaThermoData());
		 System.out.println("ThermoData is \n" +  chemgraph.getThermoData().toString());
        //int K = chemgraph.getKekule();
        int symm = chemgraph.getSymmetryNumber();
        //System.out.println("number of Kekule structures = "+K);
        System.out.println(" symmetry number = "+symm);

         Temperature T = new Temperature(298.0,"K");

         String chemicalFormula = chemgraph.getChemicalFormula();

         System.out.println(chemicalFormula + "  H=" + chemgraph.calculateH(T));
		 System.out.println();
                 ++counter;
		}

    //      Species species = Species.make(chemicalFormula, chemgraph);
          // this is equal to  System.out.println(species.toString());
    //      System.out.println(species.toStringWithoutH());
  //        species.generateResonanceIsomers();
  //        Iterator rs = species.getResonanceIsomers();
  //        while (rs.hasNext()){
  //          ChemGraph cg = (ChemGraph)rs.next();
  //          Species s = cg.getSpecies();
  //          String string = s.toStringWithoutH();
  //          System.out.print(string);
  //        }

   //       Species species = Species.make(chemicalFormula, chemgraph);
   //       Iterator iterator = species.getResonanceIsomers();
   //       System.out.println(iterator);


System.out.println("Done!\n");

};









}

