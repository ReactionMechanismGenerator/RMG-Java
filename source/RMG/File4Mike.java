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
import java.util.LinkedList;
import java.util.StringTokenizer;

import jing.chem.ChemGraph;
import jing.chem.ForbiddenStructureException;
import jing.chem.InvalidChemGraphException;
import jing.chem.Species;
import jing.chemParser.ChemParser;
import jing.chemUtil.Graph;
import jing.param.Global;
import jing.param.Temperature;
import jing.rxn.Kinetics;
import jing.rxn.Reaction;
import jing.rxn.Structure;
import jing.rxn.TemplateReactionGenerator;
import jing.rxnSys.ReactionModelGenerator;

public class File4Mike {

	public static void main(String[] args) throws InvalidChemGraphException, ForbiddenStructureException {
		initializeSystemProperties();
		// Set Global.lowTemp and Global.highTemp
		//	The values of the low/highTemp are not used in the function
		//		(to the best of my knowledge).
		//	They are necessary for the instances of additionalKinetics, 
		//	e.g. H2C*-CH2-CH2-CH3 -> H3C-CH2-*CH-CH3
		Global.lowTemperature = new Temperature(300,"K");
		Global.highTemperature = new Temperature(1500,"K");
		
		Temperature systemTemp_K = new Temperature(1200.0,"K");
		// Creating a new ReactionModelGenerator so I can set the variable temp4BestKinetics
		ReactionModelGenerator rmg = new ReactionModelGenerator();
		rmg.setTemp4BestKinetics(systemTemp_K);
		TemplateReactionGenerator rtLibrary = new TemplateReactionGenerator();
		
		int rxnCount = 0;
		// Define string to write to file
		String fileContents = "";
		
		// Open and read the input file
		try {
			FileReader fr_input = new FileReader(args[0]);
			BufferedReader br_input = new BufferedReader(fr_input);

			String line = ChemParser.readMeaningfulLine(br_input);
          
			while (line != null) {
				++rxnCount;
				System.out.println(rxnCount);
				fileContents += rxnCount + "\t";
				
				// Define variable 'speciesSet' to store the species contained in the input file
				LinkedList reactantSet = new LinkedList();
				LinkedHashSet reactants = new LinkedHashSet();
				LinkedList productSet = new LinkedList();
				// Define variable 'reactions' to store all possible rxns between the species in speciesSet
				LinkedHashSet reactions = new LinkedHashSet();
				
				// Split the string into reactants and products
	            String[] reactsANDprods = line.split("\\-->");
	            
	            // Split the reactants string into individual reactants
	            String[] reacts = reactsANDprods[0].split("[+]");
	            
	            // Make the reactant species
	            Species reactant1 = null;
	            Species reactant2 = null;
	            Species reactant3 = null;
	            
	            reactant1 = InChI2Species(reacts[0].trim());
	            reactantSet.add(reactant1);
	            reactants.add(reactant1);
	            if (reacts.length > 1) {
	            	reactant2 = InChI2Species(reacts[1].trim());
	            	reactantSet.add(reactant2);
	            	reactants.add(reactant2);
	            	if (reacts.length > 2) {
	            		reactant3 = InChI2Species(reacts[2].trim());
	            		reactantSet.add(reactant3);
	            		reactants.add(reactant3);
	            	}
	            }

	            // Split the products string into individual products
	            String[] prods = reactsANDprods[1].split("[+]");
	            
	            // Make the product species
	            Species product1 = null;
	            Species product2 = null;
	            Species product3 = null;
	            
	            product1 = InChI2Species(prods[0].trim());
	            productSet.add(product1);
	            if (prods.length > 1) {
	            	product2 = InChI2Species(prods[1].trim());
	            	productSet.add(product2);
	            	if (prods.length > 2) {
	            		product3 = InChI2Species(prods[2].trim());
	            		productSet.add(product3);
	            	}
	            }
	            
	            Structure s = new Structure(reactantSet,productSet,1);
	            
				// React the reactants
				reactions = rtLibrary.react(reactants);
				
				// Determine if reaction of interest exists in set of reactions,
				//	using the structure s as the unique identifier
				Iterator iter_rxns = reactions.iterator();
	        	while (iter_rxns.hasNext()){
	        		Reaction r = (Reaction)iter_rxns.next();
	        		
	        		if (r.getStructure().equals(s)) {
	        			fileContents += writeOutputString(r,rtLibrary);
	        			break;
	        		} else if (r.getReverseReaction().getStructure().equals(s)) {
	        			fileContents += writeOutputString(r.getReverseReaction(),rtLibrary);
	        			break;
	        		}
				}
	            
	        	fileContents += "\n";
	            line = ChemParser.readMeaningfulLine(br_input);
			}
			
        	// Write the output files
        	try{
        		File data = new File("RMG_data.txt");
        		FileWriter fw_data = new FileWriter(data);
        		fw_data.write(fileContents);
        		fw_data.close();
        	}
        	catch (IOException e) {
            	System.out.println("Could not write PopRxnsOutput.txt files");
            	System.exit(0);
            }
        	
		} catch (FileNotFoundException e) {
			System.err.println("File was not found!\n");
		} catch(IOException e) {
			System.err.println("Something wrong with ChemParser.readChemGraph");
		}
	}
	
	private static Species InChI2Species(String string) throws InvalidChemGraphException, ForbiddenStructureException {
		String InChI = Species.inchi2AdjList(string);
		Graph graph = ChemParser.readAdjList(InChI);
		return Species.make(string, ChemGraph.make(graph));
	}

	public static void initializeSystemProperties() {
		File GATPFit = new File("GATPFit");
		ChemParser.deleteDir(GATPFit);
		GATPFit.mkdir();

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
			   + "\t" + rxn_k.getComment() + "\n";//***?
		return output;
	}
	
	public static String updateListOfReactions(Kinetics rxn_k, String reverseRxnName) {
		String output = rxn_k.getAValue() + "\t" + rxn_k.getNValue()
			   + "\t" + rxn_k.getEValue() + "\t" + reverseRxnName + ": "
			   + rxn_k.getSource() + "\t" + rxn_k.getComment() + "\n";//***?
		return output;
	}
	
	public static String writeOutputString(Reaction r, TemplateReactionGenerator rtLibrary) {
		String listOfReactions = "";

		if (r.isForward()) {
			Kinetics[] allKinetics = r.getKinetics();
			for (int numKinetics=0; numKinetics<allKinetics.length; ++numKinetics) {
				listOfReactions += r.toString() + "\t" + updateListOfReactions(allKinetics[numKinetics]);
				if (allKinetics.length != 1) listOfReactions += "\tDUP\n";
			}
		}
		else if (r.isBackward()) {
			LinkedHashSet reverseReactions = new LinkedHashSet();
			Iterator iter2 = r.getStructure().getProducts();
			Species species1 = (Species)iter2.next();
			Species species2 = null;
			while (iter2.hasNext()) 
				species2 = (Species)iter2.next();
			String rxnFamilyName = r.getReverseReaction().getReactionTemplate().getName();
			reverseReactions = rtLibrary.reactSpecificFamily(species1, species2, rxnFamilyName);
			for (Iterator iter3 = reverseReactions.iterator(); iter3.hasNext();) {
				Reaction currentRxn = (Reaction)iter3.next();
				Kinetics[] allKinetics = currentRxn.getKinetics();
				if (currentRxn.getStructure() == r.getReverseReaction().getStructure()) {
					for (int numKinetics=0; numKinetics<allKinetics.length; ++numKinetics) {
						listOfReactions += currentRxn.toString() + "\t" + updateListOfReactions(allKinetics[numKinetics]);
						if (allKinetics.length != 1) listOfReactions += "\tDUP\n";
					}
				} 
			}
		}
		else {
			Kinetics[] allKinetics = r.getKinetics();
			for (int numKinetics=0; numKinetics<allKinetics.length; ++numKinetics) {
				listOfReactions += r.toString() + "\t" + updateListOfReactions(allKinetics[numKinetics]);
				if (allKinetics.length != 1) listOfReactions += "\tDUP\n";
			}
		}
		
//		if (r.hasAdditionalKinetics()) {
//			HashSet indiv_rxn = r.getAllKinetics();
//			for (Iterator iter = indiv_rxn.iterator(); iter.hasNext();) {
//				Kinetics k_rxn = (Kinetics)iter.next();
//				if (r.isForward())	listOfReactions += r.toString() + "\t" + updateListOfReactions(k_rxn) + "\tDUP";
//				else if (r.isBackward()) {
//					LinkedHashSet reverseReactions = new LinkedHashSet();
//					Iterator iter2 = r.getStructure().getProducts();
//					Species species1 = (Species)iter.next();
//					Species species2 = null;
//					while (iter2.hasNext()) 
//						species2 = (Species)iter2.next();
//					String rxnFamilyName = r.getReverseReaction().getReactionTemplate().getName();
//					reverseReactions = rtLibrary.reactSpecificFamily(species1, species2, rxnFamilyName);
//					for (Iterator iter3 = reverseReactions.iterator(); iter3.hasNext();) {
//						Reaction currentRxn = (Reaction)iter3.next();
//						if (currentRxn.getStructure() == r.getReverseReaction().getStructure()) {
//							listOfReactions += currentRxn.getReverseReaction().toString() + "\t" + updateListOfReactions(currentRxn.getReverseReaction().getFittedReverseKinetics()) + "\tDUP\tFitted Reversed Kinetics!!!";
//							break;
//						} 
//					}
//				}
//				else listOfReactions += r.toString() + "\t" + updateListOfReactions(k_rxn) + "\tDUP";
//			}
//		} else {
//			if (r.isForward()) listOfReactions += r.toString() + "\t" + updateListOfReactions(r.getKinetics());
//			else if (r.isBackward()) {
//				LinkedHashSet reverseReactions = new LinkedHashSet();
//				Iterator iter = r.getStructure().getProducts();
//				Species species1 = (Species)iter.next();
//				Species species2 = null;
//				while (iter.hasNext()) 
//					species2 = (Species)iter.next();
//				String rxnFamilyName = r.getReverseReaction().getReactionTemplate().getName();
//				reverseReactions = rtLibrary.reactSpecificFamily(species1, species2, rxnFamilyName);
//				for (Iterator iter2 = reverseReactions.iterator(); iter2.hasNext();) {
//					Reaction currentRxn = (Reaction)iter2.next();
//					if (currentRxn.getStructure() == r.getReverseReaction().getStructure()) {
//						listOfReactions += currentRxn.getReverseReaction().toString() + "\t" + updateListOfReactions(currentRxn.getReverseReaction().getFittedReverseKinetics()) + "\tFitted Reversed Kinetics!!!";
//						break;
//					}
//				}
//			}
//			else listOfReactions += r.toString() + "\t" + updateListOfReactions(r.getKinetics());
//		}
		return listOfReactions;
	}

}
