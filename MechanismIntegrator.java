import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.StringTokenizer;

import jing.chem.Species;
import jing.chemParser.ChemParser;
import jing.rxnSys.Chemkin;
import jing.rxnSys.CoreEdgeReactionModel;
import jing.rxnSys.InvalidSymbolException;
import jing.rxnSys.ReactionModelGenerator;
import jing.rxnSys.ReactionSystem;
import jing.rxnSys.ReactionTime;


public class MechanismIntegrator {
	public static void main(String[] args){
		ReactionModelGenerator rmg = new ReactionModelGenerator();
		initializeSystemProperties(args[0]);
		try {
			rmg.initializeReactionSystemForIntegration();
		} catch (InvalidSymbolException e) {
			System.err.println(e.getMessage());
        	System.exit(0);
		} catch (IOException e) {
			System.err.println(e.getMessage());
        	System.exit(0);
		}
		
		rmg.parseAllSpecies();
		rmg.parseCoreSpecies();
		rmg.parseCoreReactions();
		
		ReactionSystem reactionSystem = rmg.getReactionSystem();
		
		if (reactionSystem.getPrimaryReactionLibrary() != null){
			((CoreEdgeReactionModel)reactionSystem.getReactionModel()).addReactedSpeciesSet(reactionSystem.getPrimaryReactionLibrary().getSpeciesSet());								
			((CoreEdgeReactionModel)reactionSystem.getReactionModel()).addPrimaryReactionSet(reactionSystem.getPrimaryReactionLibrary().getReactionSet());

		}
		
		ReactionTime init = rmg.getReactionSystem().getInitialReactionTime();
		ReactionTime delt = rmg.getTimeStep();
		ReactionTime begin = init;
		ReactionTime end = begin.add(delt);
		boolean terminated = reactionSystem.isReactionTerminated();
		
		while (!terminated){
			reactionSystem.solveReactionSystem(begin, end, true, false, true);
			terminated = reactionSystem.isReactionTerminated();
			end = end.add(delt);
		}
		
		LinkedList speList = new LinkedList(reactionSystem.getReactionModel().getSpeciesSet());
		String finalOutput = "";
		finalOutput = finalOutput +"\n"+ reactionSystem.returnMoleFractionProfile(speList)+"\n";
		Chemkin.writeChemkinInputFile(reactionSystem.getReactionModel(),reactionSystem.getPresentStatus());
		
		System.out.println(finalOutput);
	}
	
	public static void initializeSystemProperties(String inputFile) {
		System.setProperty("jing.rxnSys.ReactionModelGenerator.conditionFile",inputFile);
		 String workingDir = System.getenv("RMG");
	     System.setProperty("RMG.workingDirectory", workingDir);
	     try {//svp
             String initialConditionFile = System.getProperty("jing.rxnSys.ReactionModelGenerator.conditionFile");
             if (initialConditionFile == null) {
                     System.out.println("undefined system property: jing.rxnSys.ReactionModelGenerator.conditionFile");
                     System.exit(0);
             }
             FileReader in = new FileReader(initialConditionFile);
             BufferedReader reader = new BufferedReader(in);
             String line = ChemParser.readMeaningfulLine(reader);
			 	//line = ChemParser.readMeaningfulLine(reader);
			 
             if (line.startsWith("Database")){
               StringTokenizer st = new StringTokenizer(line);
               String next = st.nextToken();
               String name = st.nextToken().trim();
               System.setProperty("jing.chem.ChemGraph.forbiddenStructureFile", workingDir + "/databases/"+name+"/forbiddenStructure/ForbiddenStructure.txt");
               System.setProperty("jing.chem.ThermoGAGroupLibrary.pathName", workingDir + "/databases/" + name+"/thermo");
               System.setProperty("jing.rxn.ReactionTemplateLibrary.pathName", workingDir + "/databases/" + name+"/kinetics/kinetics");
             }
             line = ChemParser.readMeaningfulLine(reader);
             if (line.startsWith("PrimaryThermoLibrary")){
               StringTokenizer st = new StringTokenizer(line);
               String next = st.nextToken();
               String name = st.nextToken().trim();
               String thermoDirectory = System.getProperty("jing.chem.ThermoGAGroupLibrary.pathName");
               System.setProperty("jing.chem.PrimaryThermoLibrary.pathName", thermoDirectory+"/"+name);
             }
	     }
	     catch (IOException e) {
	    	 System.err.println("Error in read in reaction system initialization file!");
	     }
	    // String name = "RMG_database";
		 //System.setProperty("jing.chem.ChemGraph.forbiddenStructureFile", workingDir + "/databases/"+name+"/forbiddenStructure/ForbiddenStructure.txt");
	     //System.setProperty("jing.chem.ThermoGAGroupLibrary.pathName", workingDir + "/databases/" + name+"/thermo");
	     //System.setProperty("jing.rxn.ReactionTemplateLibrary.pathName", workingDir + "/databases/" + name+"/kinetics/kinetics");
	     //String thermoDirectory = System.getProperty("jing.chem.ThermoGAGroupLibrary.pathName");
		//System.setProperty("jing.chem.PrimaryThermoLibrary.pathName", thermoDirectory+"/"+"primaryThermoLibrary");	   
	}
	
}

