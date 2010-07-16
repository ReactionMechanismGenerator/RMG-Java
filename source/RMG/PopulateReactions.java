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
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.StringTokenizer;

import jing.chem.ChemGraph;
import jing.chem.Species;
import jing.chemParser.ChemParser;
import jing.param.Temperature;
import jing.rxn.ArrheniusKinetics;
import jing.rxn.ArrheniusEPKinetics;
import jing.rxn.Kinetics;
import jing.rxn.LibraryReactionGenerator;
import jing.rxn.PDepIsomer;
import jing.rxn.PDepKineticsEstimator;
import jing.rxn.PDepNetwork;
import jing.rxn.PDepReaction;
import jing.rxn.Reaction;
import jing.rxn.ReactionLibrary;
import jing.rxn.TemplateReactionGenerator;
import jing.rxnSys.ConstantTM;
import jing.rxnSys.CoreEdgeReactionModel;
import jing.rxnSys.FinishController;
import jing.rxnSys.InitialStatus;
import jing.rxnSys.PressureModel;
import jing.rxnSys.RateBasedPDepRME;
import jing.rxnSys.RateBasedRME;
import jing.rxnSys.ReactionModelGenerator;
import jing.rxnSys.ReactionSystem;
import jing.rxnSys.TemperatureModel;

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
	 * 
	 *	UPDATE by MRH on 10/Jan/2010 - Fixed 2 bugs:
	 *		(1) Module did not report the correct A (pre-exponential factor) for reverse
	 *			reactions (e.g. H+CH4=CH3+H2)
	 *		(2) Module reported the same reaction multiple times (no dupliate catch)
	 *
	 *	UPDATE by MRH on 7/Apr/2010 - Restructured input file + Added pdep options
	 *		(1) In addition to reporting the non p-dep reaction rates, user will have
	 *			the option to ask for (all) the p-dep networks.
	 *		(2) The structure of the input file closely resembles the input file for
	 *			the RMG module.
	 */
	public static void main(String[] args) {
		initializeSystemProperties();
		try {
			ChemGraph.readForbiddenStructure();
		} catch (IOException e1) {
			System.err.println("PopulateReactions cannot locate forbiddenStructures.txt file");
			e1.printStackTrace();
		}
		
		// Creating a new ReactionModelGenerator so I can set the variable temp4BestKinetics
		//	and call the new readAndMakePTL and readAndMakePRL methods
		ReactionModelGenerator rmg = new ReactionModelGenerator();
		rmg.setSpeciesSeed(new LinkedHashSet());
		
		// Set Global.lowTemp and Global.highTemp
		//	The values of the low/highTemp are not used in the function
		//		(to the best of my knowledge).
		//	They are necessary for the instances of additionalKinetics, 
		//	e.g. H2C*-CH2-CH2-CH3 -> H3C-CH2-*CH-CH3
		/*
		 * 7Apr2010:
		 * The input file will now ask the user for a TemperatureModel
		 * 	and PressureModel (same as the RMG module).  The Global
		 * 	.lowTemperature and .highTemperature will automatically
		 * 	be determined
		 */
//		Global.lowTemperature = new Temperature(300,"K");
//		Global.highTemperature = new Temperature(1500,"K");

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
			
			/*
			 * Read max atom types (if they exist)
			 */
			line = rmg.readMaxAtomTypes(line,br_input);
			
			/*
			 * Read primary thermo libraries (if they exist)
			 */
            if (line.toLowerCase().startsWith("primarythermolibrary")) {
            	rmg.readAndMakePTL(br_input);
            }
            else {
            	System.err.println("PopulateReactions: Could not locate the PrimaryThermoLibrary field.\n" +
            			"Line read was: " + line);
            	System.exit(0);
            }
            line = ChemParser.readMeaningfulLine(br_input);
            
            // Read primary transport library
            if (line.toLowerCase().startsWith("primarytransportlibrary"))
            	rmg.readAndMakePTransL(br_input);
            else {
            	System.err.println("PopulateReactions: Could not locate the PrimaryTransportLibrary field.\n" +
            			"Line read was: " + line);
            	System.exit(0);
            }
			
            /*
             * Read the temperature model (must be of length one)
             */
            line = ChemParser.readMeaningfulLine(br_input);
			rmg.createTModel(line);
//			Temperature systemTemp = null;
//			if (!line.startsWith("Temperature")) {
//				System.err.println("Error reading input file: Could not locate System Temperature.\n" +
//						"The first line of the input file should read: \"Temperature: Value (Units)\"");
//			} else {
//				StringTokenizer st = new StringTokenizer(line);
//				String dummyString = st.nextToken();	// This token should be "Temperature:"
//				systemTemp = new Temperature(Double.parseDouble(st.nextToken()),ChemParser.removeBrace(st.nextToken()));
//			}
//			Temperature systemTemp_K = new Temperature(systemTemp.getK(),"K");			
//			listOfReactions += "System Temperature: " + systemTemp_K.getK() + "K\n\n";
			if (rmg.getTempList().size() > 1) {
				System.out.println("Please list only one temperature in the TemperatureModel field.");
				System.exit(0);
			}
			
			/*
			 * Read the pressure model (must be of length 1)
			 */
			line = ChemParser.readMeaningfulLine(br_input);
			rmg.createPModel(line);
			if (rmg.getPressList().size() > 1) {
				System.out.println("Please list only one pressure in the PressureModel field.");
				System.exit(0);
			}
			
			/*
			 * Read the solvation field (if present)
			 */
			line = ChemParser.readMeaningfulLine(br_input);
            StringTokenizer st = new StringTokenizer(line);
            // The first line should start with "Solvation", otherwise do nothing and display a message to the user
            if (st.nextToken().startsWith("Solvation")) {
            	line = st.nextToken().toLowerCase();
            	// The options for the "Solvation" field are "on" or "off" (as of 18May2009), otherwise do nothing and display a message to the user
            	// Note: I use "Species.useInChI" because the "Species.useSolvation" updates were not yet committed.

            	if (line.equals("on")) {
            		Species.useSolvation = true;
            		//                  rmg.setUseDiffusion(true);
            		listOfReactions += "Solution-phase chemistry!\n\n";
            	} else if (line.equals("off")) {
            		Species.useSolvation = false;
            		//                  rmg.setUseDiffusion(false);
            		listOfReactions += "Gas-phase chemistry.\n\n";
            	} else {
            		System.out.println("Error in reading input.txt file:\nThe field 'Solvation' has the options 'on' or 'off'." +
            				"\nPopulateReactions does not recognize: " + line);
            		return;
            	}
            	line = ChemParser.readMeaningfulLine(br_input);
            }

            /*
             * Read in the species (name, concentration, adjacency list)
             */
            if (line.toLowerCase().startsWith("speciesstatus")) {
            	LinkedHashMap lhm = new LinkedHashMap();
            	lhm = rmg.populateInitialStatusListWithReactiveSpecies(br_input);
            	speciesSet.addAll(lhm.values());
            }
            
            /*
             * Read in the inert gas (name, concentration)
             */
            line = ChemParser.readMeaningfulLine(br_input);
            if (line.toLowerCase().startsWith("bathgas")) {
            	rmg.populateInitialStatusListWithInertSpecies(br_input);
            }
            
            /*
             * Read in the p-dep options
             */
            line = ChemParser.readMeaningfulLine(br_input);
            if (line.toLowerCase().startsWith("spectroscopicdata")) {
            	rmg.setSpectroscopicDataMode(line);
            	line = ChemParser.readMeaningfulLine(br_input);
            	line = rmg.setPressureDependenceOptions(line,br_input);
            }
            
            /*
             * Read primary kinetic libraries (if they exist)
             */
            if (line.toLowerCase().startsWith("primarykineticlibrary")) {
            	rmg.readAndMakePKL(br_input);
            }
                        
           else {
              	System.err.println("PopulateReactions: Could not locate the PrimaryKineticLibrary field." +
             			"Line read was: " + line);
             	System.exit(0);
             }
            
            line = ChemParser.readMeaningfulLine(br_input);
            if (line.toLowerCase().startsWith("reactionlibrary")) {
            	rmg.readAndMakeReactionLibrary(br_input);
            }
            else {
              	System.err.println("PopulateReactions: Could not locate the ReactionLibrary field." +
             			"Line read was: " + line);
             	System.exit(0);
             }
            
            
            /*
             * MRH 2Apr2010:
             * Allow user the option to print "verbose" comments
             */
            /*
             * Read in verbosity field (if it exists)
             */
            line = ChemParser.readMeaningfulLine(br_input);
            if (line != null && line.toLowerCase().startsWith("verbose")) {
            	StringTokenizer st2 = new StringTokenizer(line);
            	String tempString = st2.nextToken();
            	tempString = st2.nextToken();
            	tempString = tempString.toLowerCase();
            	if (tempString.equals("on") || tempString.equals("true") || tempString.equals("yes"))
            		ArrheniusKinetics.setVerbose(true);
            }
            
            /*
             * MRH 27Feb:
             * Allowing PopulateReactions module to read/use Primary Kinetic/Thermo Libraries
             * 
             * The user does not have to specify any primary kinetic / thermo libraries.
             * 	However, like the input file to RMG, the fields must still be present
             */
//            if (line.toLowerCase().startsWith("primarythermolibrary")) {
//            	rmg.readAndMakePTL(br_input);
//            }
//            else {
//            	System.err.println("PopulateReactions: Could not locate the PrimaryThermoLibrary field." +
//            			"Line read was: " + line);
//            	System.exit(0);
//            }
            
//            line = ChemParser.readMeaningfulLine(br_input);
//            if (line.toLowerCase().startsWith("primarykineticlibrary")) {
//            	rmg.readAndMakePKL(br_input);
//            }
//            else {
//            	System.err.println("PopulateReactions: Could not locate the PrimaryKineticLibrary field." +
//            			"Line read was: " + line);
//            	System.exit(0);
//            }
            
//            line = ChemParser.readMeaningfulLine(br_input);
//            BathGas bathgas = null;
//            String inertgasname = "";
//            if (line.toLowerCase().startsWith("spectroscopicdata")) {
//            	rmg.setSpectroscopicDataMode(line);
//            	line = ChemParser.readMeaningfulLine(br_input);
//            	line = rmg.setPressureDependenceOptions(line,br_input);
//            	if (PDepNetwork.generateNetworks) {
//	            	StringTokenizer bathGas_st = new StringTokenizer(line);
//	            	String name = bathGas_st.nextToken();
//	            	name = bathGas_st.nextToken();
//	            	inertgasname = name;
//	            	bathgas = new BathGas(name);
//	            	line = ChemParser.readMeaningfulLine(br_input);
//            	}
//            }
            
//            LinkedHashMap dummyLHM = new LinkedHashMap();
//            
//            // Read in all of the species in the input file
//			while (line != null) {
//				// The first line of a new species is the user-defined name
//				String speciesName = line;
//				// The remaining lines are the graph
//				Graph g = ChemParser.readChemGraph(br_input);
//				// Make the ChemGraph, assuming it does not contain a forbidden structure
//				ChemGraph cg = null;
//				try {
//					cg = ChemGraph.make(g);
//				} catch (ForbiddenStructureException e) {
//					System.out.println("Error in reading graph: Graph contains a forbidden structure.\n" + g.toString());
//					System.exit(0);
//				}
//				// Make the species
//				Species species = Species.make(speciesName,cg);
//				// Add the new species to the set of species
//				//int speciesCount = speciesSet.size();
//				speciesSet.add(species);
//    			SpeciesStatus ss = new SpeciesStatus(species,1,1e-12,0);
//    			dummyLHM.put(species, ss);
//				//if (speciesSet.size() != speciesCount+1){
//				//	System.out.println(speciesName);
//				//	Iterator iter=speciesSet.iterator();
//				//	while (iter.hasNext()) {
//				//		Species spcs = (Species)iter.next();
//				//		if (spcs.equals(species))
//				//			System.out.println(spcs.toString());
//				//	}
//				//}
//
//				// Read the next line of the input file
//				line = ChemParser.readMeaningfulLine(br_input);
//			}
//			InitialStatus is = new InitialStatus(dummyLHM,null,null);
//			is.putInertGas(inertgasname,1e-6);
			
			// Set the user's input temperature
            LinkedList tempList = rmg.getTempList();
            Temperature systemTemp = ((ConstantTM)tempList.get(0)).getTemperature();
			rmg.setTemp4BestKinetics(systemTemp);
			TemplateReactionGenerator rtLibrary = new TemplateReactionGenerator();
			ReactionLibrary  RL = rmg.getReactionLibrary();
			LibraryReactionGenerator lrg1 =new LibraryReactionGenerator(RL);
			// React all species with each other
			//reactions = rtLibrary.react(speciesSet);
			
			reactions = lrg1.react(speciesSet);
			System.out.println("Reactions From LRG"+reactions);
			reactions.addAll(rtLibrary.react(speciesSet));
			/*
			 * MRH 28Feb
			 * The next line will use the new PrimaryKineticLibrary
			 */
			//reaction.addAll(getLibraryReactionGenerator().react(speciesSet));
			if (!(rmg.getReactionModelEnlarger() instanceof RateBasedRME))	{
				CoreEdgeReactionModel cerm = new CoreEdgeReactionModel(speciesSet, reactions);
				rmg.setReactionModel(cerm);
				rmg.setReactionGenerator(rtLibrary);
				
				ReactionSystem rs = new ReactionSystem((TemperatureModel)rmg.getTempList().get(0),
						(PressureModel)rmg.getPressList().get(0),
						rmg.getReactionModelEnlarger(),
						new FinishController(),
						null,
						rmg.getPrimaryKineticLibrary(),
						rmg.getReactionGenerator(),
						speciesSet,
						(InitialStatus)rmg.getInitialStatusList().get(0),
						rmg.getReactionModel(),
						rmg.getLibraryReactionGenerator(),
						0,
						"GasPhase");
				
				PDepNetwork.reactionModel = rmg.getReactionModel();
				PDepNetwork.reactionSystem = rs;
				
				// If the reaction structure is A + B = C + D, we are not concerned w/pdep
				Iterator iter = reactions.iterator();
				LinkedHashSet nonPdepReactions = new LinkedHashSet();
	        	while (iter.hasNext()){
	        		Reaction r = (Reaction)iter.next();
	        		if (r.getReactantNumber() < 2 || r.getProductNumber() < 2){
						cerm.categorizeReaction(r.getStructure());
						PDepNetwork.addReactionToNetworks(r);
					}
	        		else {
	        			nonPdepReactions.add(r);
	        		}
				}
	        	
	        	// Run fame calculation
	        	PDepKineticsEstimator pDepKineticsEstimator = 
					((RateBasedPDepRME) rmg.getReactionModelEnlarger()).getPDepKineticsEstimator();
	        	
	        	for (int numNetworks=0; numNetworks<PDepNetwork.getNetworks().size(); ++numNetworks) {
	        		PDepNetwork pdepnetwork = PDepNetwork.getNetworks().get(numNetworks);
		        	pDepKineticsEstimator.runPDepCalculation(pdepnetwork,rs,cerm);
		        	
		        	LinkedList<PDepReaction> indivPDepRxns = pdepnetwork.getNetReactions();
		        	for (int numPDepRxns=0; numPDepRxns<indivPDepRxns.size(); numPDepRxns++) {
		        		listOfReactions += indivPDepRxns.get(numPDepRxns).toRestartString(systemTemp);
		        	}
		        	LinkedList<PDepReaction> nonIncludedRxns = pdepnetwork.getNonincludedReactions();
		        	for (int numNonRxns=0; numNonRxns<nonIncludedRxns.size(); ++numNonRxns) {
		        		listOfReactions += nonIncludedRxns.get(numNonRxns).toRestartString(systemTemp);
		        	}
		        	LinkedList<PDepIsomer> allpdepisomers = pdepnetwork.getIsomers();
		        	for (int numIsomers=0; numIsomers<allpdepisomers.size(); ++numIsomers) {
		        		LinkedList species = allpdepisomers.get(numIsomers).getSpeciesList();
		        		speciesSet.addAll(species);		        		
		        	}
	        	}
	        	reactions = nonPdepReactions;
			}			
			
			// Some of the reactions may be duplicates of one another 
			//	(e.g. H+CH4=CH3+H2 as a forward reaction and reverse reaction)
			//	Create new LinkedHashSet which will store the non-duplicate rxns
			LinkedHashSet nonDuplicateRxns = new LinkedHashSet();
			int Counter = 0;
			
			Iterator iter_rxns = reactions.iterator();
        	while (iter_rxns.hasNext()){
        		++Counter;
        		Reaction r = (Reaction)iter_rxns.next();
        		
        		// The first reaction is not a duplicate of any previous reaction
        		if (Counter == 1) {
        			nonDuplicateRxns.add(r);
        			listOfReactions += writeOutputString(r,rtLibrary);
            		speciesSet.addAll(r.getProductList());
        		}
        		
        		// Check whether the current reaction (or its reverse) has the same structure
        		//	of any reactions already reported in the output
        		else {
        			Iterator iterOverNonDup = nonDuplicateRxns.iterator();
        			boolean dupRxn = false;
	        		while (iterOverNonDup.hasNext()) {
	        			Reaction temp_Reaction = (Reaction)iterOverNonDup.next();
	        			if (r.getStructure() == temp_Reaction.getStructure()) {
	        				dupRxn = true;
	        				break;
	        			} else if (r.hasReverseReaction()) {
	        				if (r.getReverseReaction().getStructure() == temp_Reaction.getStructure()) {
	        					dupRxn = true;
	        					break;
	        				}
	        			}
	        		}
	        		if (!dupRxn) {
        				nonDuplicateRxns.add(r);
        				// If Reaction is Not a Library Reaction
        				listOfReactions += writeOutputString(r,rtLibrary);
                		speciesSet.addAll(r.getProductList());
	        		}
        		}
			}
        	
        	Iterator iter_species = speciesSet.iterator();
        	// Define dummy integer 'i' so our getChemGraph().toString()
        	//	call only returns the graph
        	int i = 0;
        	while (iter_species.hasNext()) {
        		Species species = (Species)iter_species.next();
        		listOfSpecies += species.getName()+"("+species.getID()+")\n" +
        			species.getChemGraph().toStringWithoutH(i) + "\n";
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
		RMG.globalInitializeSystemProperties();
		
		File GATPFit = new File("GATPFit");
		ChemParser.deleteDir(GATPFit);
		GATPFit.mkdir();
		File frankie = new File("frankie");
		ChemParser.deleteDir(frankie);
		frankie.mkdir();
		File fame = new File("fame");
		ChemParser.deleteDir(fame);
		fame.mkdir();
	};
	
	public static String updateListOfReactions(Kinetics rxn_k, double H_rxn) {
		double Ea = 0.0;
		if (rxn_k instanceof ArrheniusEPKinetics){
		    Ea = ((ArrheniusEPKinetics)rxn_k).getEaValue(H_rxn);
		}
		else{
		    Ea = rxn_k.getEValue();
		}
		String output = rxn_k.getAValue() + "\t" + rxn_k.getNValue()
			   + "\t" + Ea + "\t" + rxn_k.getSource()
			   + "\t" + rxn_k.getComment() 
			   + "\tdeltaHrxn(T=298K) = " + H_rxn + " kcal/mol\n";
		return output;
	}
	
	public static String writeOutputString(Reaction r, TemplateReactionGenerator rtLibrary) {
		String listOfReactions = "";

		Temperature stdtemp = new Temperature(298,"K");
		double Hrxn = r.calculateHrxn(stdtemp);
		
		// If r Reaction is from Reaction Library add it to list of reaction append its kinetics and return
        String source = r.getKineticsSource(0);
        if (source == null){
    		// If source is null I am assuming that its not a Reaction from Reaction Library or Seed Mechanism
    		source = "TemplateReaction:";
    	}
        StringTokenizer st = new StringTokenizer(source,":");
    	String reaction_type = st.nextToken();
    	if (reaction_type.equals("ReactionLibrary") ){
    		Kinetics[] allKinetics = r.getKinetics();
			for (int numKinetics=0; numKinetics<allKinetics.length; ++numKinetics) {
				listOfReactions += r.toString() + "\t" + updateListOfReactions(allKinetics[numKinetics], Hrxn);
				if (allKinetics.length != 1) listOfReactions += "\tDUP\n";
			}
		return 	listOfReactions;
       	}
        
		if (r.isForward()) {
			Kinetics[] allKinetics = r.getKinetics();
			for (int numKinetics=0; numKinetics<allKinetics.length; ++numKinetics) {
				listOfReactions += r.toString() + "\t" + updateListOfReactions(allKinetics[numKinetics], Hrxn);
				if (allKinetics.length != 1) listOfReactions += "\tDUP\n";
			}
		}
		else if (r.isBackward()) {
			LinkedHashSet reverseReactions = new LinkedHashSet();
			Iterator iter2 = r.getStructure().getProducts();
			Species species1 = (Species)iter2.next();
			Species species2 = species1;
			while (iter2.hasNext()) 
				species2 = (Species)iter2.next();
			String rxnFamilyName = r.getReverseReaction().getReactionTemplate().getName();
			LinkedHashSet speciesHashSet = new LinkedHashSet();
			speciesHashSet.add(species1);
			reverseReactions = rtLibrary.react(speciesHashSet, species2, rxnFamilyName);
			for (Iterator iter3 = reverseReactions.iterator(); iter3.hasNext();) {
				Reaction currentRxn = (Reaction)iter3.next();
				if (currentRxn.getStructure() == r.getReverseReaction().getStructure()) {
					Kinetics[] allKinetics = currentRxn.getKinetics();
					for (int numKinetics=0; numKinetics<allKinetics.length; ++numKinetics) {
						listOfReactions += currentRxn.toString() + "\t" + updateListOfReactions(allKinetics[numKinetics], -Hrxn);
						if (allKinetics.length != 1) listOfReactions += "\tDUP\n";
					}
				} 
			}
		}
		else {
			Kinetics[] allKinetics = r.getKinetics();
			for (int numKinetics=0; numKinetics<allKinetics.length; ++numKinetics) {
				listOfReactions += r.toString() + "\t" + updateListOfReactions(allKinetics[numKinetics], Hrxn);
				if (allKinetics.length != 1) listOfReactions += "\tDUP\n";
			}
		}
//		if (r.hasAdditionalKinetics()) {
//			HashSet indiv_rxn = r.getAllKinetics();
//			for (Iterator iter = indiv_rxn.iterator(); iter.hasNext();) {
//				Kinetics k_rxn = (Kinetics)iter.next();
//				if (r.isForward())	listOfReactions += r.toString() + "\t" + updateListOfReactions(k_rxn) + "\tDUP\n";
//				//else if (r.isBackward()) listOfReactions += r.getReverseReaction().toString() + "\t" + updateListOfReactions(k_rxn) + "\tDUP\n";
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
//						if (currentRxn.getStructure() == r.getReverseReaction().getStructure())
//							listOfReactions += currentRxn.toString() + "\t" + updateListOfReactions(currentRxn.getKinetics()) + "\tDUP\n"; 
//					}
//					//listOfReactions += r.getReverseReaction().toString() + "\t" + updateListOfReactions(r.getReverseReaction().getKinetics());
//				}
//				else listOfReactions += r.toString() + "\t" + updateListOfReactions(k_rxn) + "\tDUP\n";
//			}
//		} else {
//			if (r.isForward()) listOfReactions += r.toString() + "\t" + updateListOfReactions(r.getKinetics());
//			//else if (r.isBackward()) listOfReactions += r.toString() + "\t" + updateListOfReactions(r.getFittedReverseKinetics(),r.getReactionTemplate().getName());
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
//					if (currentRxn.getStructure() == r.getReverseReaction().getStructure())
//						listOfReactions += currentRxn.toString() + "\t" + updateListOfReactions(currentRxn.getKinetics()); 
//				}
//				//listOfReactions += r.getReverseReaction().toString() + "\t" + updateListOfReactions(r.getReverseReaction().getKinetics());
//			}
//			else listOfReactions += r.toString() + "\t" + updateListOfReactions(r.getKinetics());
//		}
		return listOfReactions;
	}


}
