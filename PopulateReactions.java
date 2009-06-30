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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.StringTokenizer;

import jing.chem.ChemGraph;
import jing.chem.ForbiddenStructureException;
import jing.chem.Species;
import jing.chemParser.ChemParser;
import jing.chemUtil.Graph;
import jing.param.Global;
import jing.param.Temperature;
import jing.rxn.Kinetics;
import jing.rxn.Reaction;
import jing.rxn.TemplateReactionGenerator;
import jing.rxnSys.ReactionModelGenerator;

public class PopulateReactions {
	/**
	 * Generates a list of all possible reactions, and their modified Arrhenius
	 * 	parameters, between all species supplied in the input file.
	 * 
	 * The input file's first line should specify the system temperature.  The
	 * 	line should be formatted as follows:
	 * 		Temperature: 500 (K)
	 * The input (.txt) file should contain a list of species, with the first
	 * 	line of each species being a user-defined name for the species and the 
	 * 	remaining lines	containing the graph (in the form of an adjacency list).
	 *  There is no limit on the number of ChemGraphs the user may supply.
	 *  
	 * The output of the function is two .txt files: PopRxnsOutput_rxns.txt and
	 * 	PopRxnsOutput_spcs.txt.  The first contains the list of reactions (including
	 * 	the modified Arrhenius parameters and RMG-generated comments) and the second
	 * 	the list of species (including the ChemGraph) that are involved in the list
	 * 	of reactions.
	 */
	public static void main(String[] args) {
		initializeSystemProperties();
		// Set Global.lowTemp and Global.highTemp
		//	The values of the low/highTemp are not used in the function
		//		(to the best of my knowledge).
		//	They are necessary for the instances of additionalKinetics, 
		//	e.g. H2C*-CH2-CH2-CH3 -> H3C-CH2-*CH-CH3
		Global.lowTemperature = new Temperature(300,"K");
		Global.highTemperature = new Temperature(1500,"K");
		// Define an integer to count the number of sets of duplicate kinetics
		int numDupKinetics = 0;
		// Define variable 'speciesSet' to store the species contained in the input file
		LinkedHashSet speciesSet = new LinkedHashSet();
		// Define variable 'reactions' to store all possible rxns between the species in speciesSet
		LinkedHashSet reactions = new LinkedHashSet();
		
		// Define two string variables 'listOfReactions' and 'listOfSpecies'
		//	These strings will hold the list of rxns (including the structure,
		//	modified Arrhenius parameters, and source/comments) and the list of
		//	species (including the chemkin name and graph), respectively
		String listOfReactions = "Arrhenius 'A' parameter has units of: mol,cm3,s\n" +
			"Arrhenius 'n' parameter is unitless\n" +
			"Arrhenius 'E' parameter has units of: kcal/mol\n\n";
		String listOfSpecies = "";
		
		// Open and read the input file
		try {
			FileReader fr_input = new FileReader(args[0]);
			BufferedReader br_input = new BufferedReader(fr_input);
			// Read in the first line of the input file
			//	This line should hold the temperature of the system, e.g.
			//		Temperature: 500 (K)
			String line = ChemParser.readMeaningfulLine(br_input);
			Temperature systemTemp = null;
			if (!line.startsWith("Temperature")) {
				System.err.println("Error reading input file: Could not locate System Temperature.\n" +
						"The first line of the input file should read: \"Temperature: Value (Units)\"");
			} else {
				StringTokenizer st = new StringTokenizer(line);
				String dummyString = st.nextToken();	// This token should be "Temperature:"
				systemTemp = new Temperature(Double.parseDouble(st.nextToken()),ChemParser.removeBrace(st.nextToken()));
			}
			Temperature systemTemp_K = new Temperature(systemTemp.getK(),"K");
			// Creating a new ReactionModelGenerator so I can set the variable temp4BestKinetics
			ReactionModelGenerator rmg = new ReactionModelGenerator();
			rmg.setTemp4BestKinetics(systemTemp_K);
			TemplateReactionGenerator rtLibrary = new TemplateReactionGenerator();
			listOfReactions += "System Temperature: " + systemTemp_K.getK() + "K\n\n";
			line = ChemParser.readMeaningfulLine(br_input);
			while (line != null) {
				// The first line of a new species is the user-defined name
				String speciesName = line;
				// The remaining lines are the graph
				Graph g = ChemParser.readChemGraph(br_input);
				// Make the ChemGraph, assuming it does not contain a forbidden structure
				ChemGraph cg = null;
				try {
					cg = ChemGraph.make(g);
				} catch (ForbiddenStructureException e) {
					System.out.println("Error in reading graph: Graph contains a forbidden structure.\n" + g.toString());
					System.exit(0);
				}
				// Make the species
				Species species = Species.make(speciesName,cg);
				// Add the new species to the set of species
				speciesSet.add(species);
				// Define a dummy hash set 'new_reactions' to hold the list
				//	of rxns for the new species against the species seed
				LinkedHashSet new_reactions = new LinkedHashSet();
				// React the new species with the set of species (including itself)
				new_reactions = rtLibrary.react(speciesSet,species);
				reactions.addAll(new_reactions);
				// Read the next line of the input file
				line = ChemParser.readMeaningfulLine(br_input);
			}
			// Iterate over all rxns, abstracting the necessary information,
			//	and store the results in the 'listOfReactions' string
			Iterator iter_rxns = reactions.iterator();
        	while (iter_rxns.hasNext()){
        		Reaction r = (Reaction)iter_rxns.next();
        		if (r.hasAdditionalKinetics()) {
        			++numDupKinetics;
        			HashSet indiv_rxn = r.getAllKinetics();
        			for (Iterator iter = indiv_rxn.iterator(); iter.hasNext();) {
        				Kinetics k_rxn = (Kinetics)iter.next();
        				if (r.isForward())	listOfReactions += r.toString() + "\t" + updateListOfReactions(k_rxn) + "\tDUP " + numDupKinetics + "\n";
        				else if (r.isBackward()) listOfReactions += r.getReverseReaction().toString() + "\t" + updateListOfReactions(k_rxn) + "\tDUP " + numDupKinetics + "\n";
        				else listOfReactions += r.toString() + "\t" + updateListOfReactions(k_rxn) + "\tDUP " + numDupKinetics + "\n";
        			}
        		} else {
        			if (r.isForward()) listOfReactions += r.toString() + "\t" + updateListOfReactions(r.getKinetics());
        			//else if (r.isBackward()) listOfReactions += r.toString() + "\t" + updateListOfReactions(r.getFittedReverseKinetics(),r.getReactionTemplate().getName());
        			else if (r.isBackward()) listOfReactions += r.getReverseReaction().toString() + "\t" + updateListOfReactions(r.getReverseReaction().getKinetics());
        			else listOfReactions += r.toString() + "\t" + updateListOfReactions(r.getKinetics());
        		}

        		// Add the products of the reactions to the list of species
        		//	hash set.  The reactants of each reaction are already
        		//	present.  This list will allow us to generate the graphs
        		//	for all species involved in rxns.
        		speciesSet.addAll(r.getProductList());
			}
        	Iterator iter_species = speciesSet.iterator();
        	// Define dummy integer 'i' so our getChemGraph().toString()
        	//	call only returns the graph
        	int i = 0;
        	while (iter_species.hasNext()) {
        		Species species = (Species)iter_species.next();
        		listOfSpecies += species.getChemkinName() + "\n" +
        			species.getChemGraph().toString(i) + "\n";
        	}
        	// Write the output files
        	try{
        		File rxns = new File("PopRxnsOutput_rxns.txt");
        		FileWriter fw_rxns = new FileWriter(rxns);
        		fw_rxns.write(listOfReactions);
        		fw_rxns.close();
        		File spcs = new File("PopRxnsOutput_spcs.txt");
        		FileWriter fw_spcs = new FileWriter(spcs);
        		fw_spcs.write(listOfSpecies);
        		fw_spcs.close();
        	}
        	catch (IOException e) {
            	System.out.println("Could not write PopRxnsOutput.txt files");
            	System.exit(0);
            }
        	// Display to the user that the program was successful and also
        	//	inform them where the results may be located
			System.out.println("Reaction population complete. Results are stored" 
					+ " in PopRxnsOutput_rxns.txt and PopRxnsOutput_spcs.txt");
		} catch (FileNotFoundException e) {
			System.err.println("File was not found!\n");
		} catch(IOException e) {
			System.err.println("Something wrong with ChemParser.readChemGraph");
		}
	}
	
	public static void initializeSystemProperties() {
		String name= "RMG_database";
		String workingDir = System.getenv("RMG");
		System.setProperty("RMG.workingDirectory", workingDir);
		System.setProperty("jing.chem.ChemGraph.forbiddenStructureFile",workingDir + "/databases/" + name + "/forbiddenStructure/forbiddenStructure.txt");
		System.setProperty("jing.chem.ThermoGAGroupLibrary.pathName",	workingDir + "/databases/" + name + "/thermo");
		System.setProperty("jing.rxn.ReactionTemplateLibrary.pathName",	workingDir + "/databases/" + name + "/kinetics");
	};
	
	public static String updateListOfReactions(Kinetics rxn_k) {
		String output = rxn_k.getAValue() + "\t" + rxn_k.getNValue()
			   + "\t" + rxn_k.getEValue() + "\t" + rxn_k.getSource() 
			   + "\t" + rxn_k.getComment() + "\n";
		return output;
	}
	
	public static String updateListOfReactions(Kinetics rxn_k, String reverseRxnName) {
		String output = rxn_k.getAValue() + "\t" + rxn_k.getNValue()
			   + "\t" + rxn_k.getEValue() + "\t" + reverseRxnName + ": "
			   + rxn_k.getSource() + "\t" + rxn_k.getComment() + "\n";
		return output;
	}

}
