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


import java.util.*;
import java.io.*;
import java.text.SimpleDateFormat;

import jing.chem.*;
import jing.chemParser.*;
import jing.param.*;
import jing.chemUtil.*;
import jing.rxn.*;
import jing.rxnSys.*;
import jing.mathTool.*;

//----------------------------------------------------------------------------
// MainRMG.java
//----------------------------------------------------------------------------

public class RMG {
	
    //## configuration RMG::RMG
    public static void main(String[] args) {
		
		// Initialize the logger
        Logger.initialize();

        // Log the RMG header
        Logger.logHeader();
        
        try {

            // Record the time at which RMG was started (so we can periodically
            // print the elapsed time)
            Global.tAtInitialization = System.currentTimeMillis();

            // Initialize some properties
            globalInitializeSystemProperties();

            // Create the working folders for various RMG components
            // Note that some of these folders are only used when their
            // corresponding features are activated in the condition file
            createFolder("chemkin", true);
            createFolder("Restart", true);
            createFolder("GATPFit", true);
            createFolder("ODESolver", true);
            createFolder("fame", true);
            createFolder("frankie", true);
            createFolder("InChI", true);
            createFolder("Pruning", true);
            createFolder("2Dmolfiles", true);   // Not sure if we should be deleting this
            createFolder("3Dmolfiles", true);   // Not sure if we should be deleting this
            createFolder("QMfiles", false);     // Preserving QM files between runs will speed things up considerably
            
            // The only parameter should be the path to the condition file
            String inputfile = args[0];
            System.setProperty("jing.rxnSys.ReactionModelGenerator.conditionFile",inputfile);
            // Echo the contents of the condition file into the log
            logInputFile(inputfile);
            
            // Read the "Database" keyword from the condition file
            // This database supercedes the default RMG_database if present
            FileReader in = new FileReader(inputfile);
            BufferedReader reader = new BufferedReader(in);
            String line = ChemParser.readMeaningfulLine(reader, true);
            if (line.startsWith("Database")) 
                extractAndSetDatabasePath(line);

            // Generate the model!
            ReactionModelGenerator rmg = new ReactionModelGenerator();
            rmg.modelGeneration();

            // Save the resulting model to Final_Model.txt
            writeFinalModel(rmg);

       }
       catch (Exception e) {
           // Any unhandled exception will land here
           // We assume that these were critical errors
           e.printStackTrace();
           Logger.critical(e.getMessage());
       }

       // Finish the logger
       Logger.finish();
       
    }

    public static void writeFinalModel(ReactionModelGenerator rmg) {
        CoreEdgeReactionModel cerm = (CoreEdgeReactionModel)rmg.getReactionModel();

        long tAtInitialization = Global.tAtInitialization;
        
		//rs.reduceModel();
        //Write the final model output in a separete Final_Model file
        String finalOutput = "";
        //finalOutput = finalOutput + "\n\\\\\\\\\\\\\\\\\\\\\\\\\\    Final Reaction Model Output    \\\\\\\\\\\\\\\\\\\\\\\\\\";
        //finalOutput = finalOutput +"\n"+ cerm.returnPDepModel(rs.getPresentStatus())+"\n";

        //finalOutput = finalOutput + "Model Edge:";
        //finalOutput = finalOutput + cerm.getEdge().getSpeciesNumber()+" ";
        //finalOutput = finalOutput + " Species; ";
        //finalOutput = finalOutput + cerm.getEdge().getReactionNumber()+" ";
        //finalOutput = finalOutput + " Reactions.";

        LinkedList speList = new LinkedList(rmg.getReactionModel().getSpeciesSet());

        //11/2/07 gmagoon: changing to loop over all reaction systems
        LinkedList rsList = rmg.getReactionSystemList();
        for (int i = 0; i < rsList.size(); i++) {
            ReactionSystem rs = (ReactionSystem) rsList.get(i);
            finalOutput = finalOutput + "Reaction system " + (i + 1) + " of " + rsList.size() + "\n";//11/4/07 gmagoon: fixing "off by one" error
            finalOutput = finalOutput + "\\\\\\\\\\\\\\\\\\\\\\\\\\\\    Concentration Profile Output    \\\\\\\\\\\\\\\\\\\\\\\\\\";
            finalOutput = finalOutput + "\n" + rs.returnConcentrationProfile(speList) + "\n";
            if (rmg.getError()) {//svp
                finalOutput = finalOutput + "\n";
                finalOutput = finalOutput + "Upper Bounds:" + "\n";
                finalOutput = finalOutput + rs.printUpperBoundConcentrations(speList) + "\n";
                finalOutput = finalOutput + "Lower Bounds:" + "\n";
                finalOutput = finalOutput + rs.printLowerBoundConcentrations(speList) + "\n";
            }
            finalOutput = finalOutput + "\\\\\\\\\\\\\\\\\\\\\\\\\\\\    Mole Fraction Profile Output    \\\\\\\\\\\\\\\\\\\\\\\\\\";
            finalOutput = finalOutput + "\n" + rs.returnMoleFractionProfile(speList) + "\n";
            finalOutput = finalOutput + rs.printOrderedReactions() + "\n";
            if (rmg.getSensitivity()) {//svp
                LinkedList importantSpecies = rmg.getSpeciesList();
                finalOutput = finalOutput + "Sensitivity Analysis:";
                finalOutput = finalOutput + rs.printSensitivityCoefficients(speList, importantSpecies) + "\n";
                finalOutput = finalOutput + "\n";
                finalOutput = finalOutput + rs.printSensitivityToThermo(speList, importantSpecies) + "\n";
                finalOutput = finalOutput + "\n";
                finalOutput = finalOutput + rs.printMostUncertainReactions(speList, importantSpecies) + "\n\n";

            }
            finalOutput = finalOutput + rs.returnReactionFlux() + "\n";
            long end = System.currentTimeMillis();
            double min = (end - tAtInitialization) / 1.0E3 / 60.0;
            finalOutput = finalOutput + "Running Time is: " + String.valueOf(min) + " minutes.";


            try {
                File finalmodel = new File("Final_Model.txt");
                FileWriter fw = new FileWriter(finalmodel);
                fw.write(finalOutput);
                fw.close();
            } catch (IOException e) {
                System.out.println("Could not write Final_Model.txt");
                System.exit(0);
            }
        }
    }

    public static void globalInitializeSystemProperties() {
    	// Set the default locale to US (since RMG is developed there)
        Locale.setDefault(Locale.US);

        // Get the value of the RMG environment variable; fail if not set!
        String workingDir = System.getenv("RMG");
        if (workingDir == null) throw new RuntimeException("The RMG environment variable is not defined.");
        System.setProperty("RMG.workingDirectory", workingDir);
		
        // Set the default database
        setDatabasePaths(workingDir + "/databases/RMG_database");

    }

	public static void setDatabasePaths(String database_path) {
		// String database_path = workingDir + "/databases/" + name
		System.setProperty("jing.chem.ChemGraph.forbiddenStructureFile",   database_path +"/ForbiddenStructures.txt");
		System.setProperty("jing.chem.ThermoGAGroupLibrary.pathName",      database_path +"/thermo_groups");
		System.setProperty("jing.chem.ThermoReferenceLibrary.pathName",    database_path +"/thermo_libraries");
		System.setProperty("jing.chem.FrequencyDatabase.pathName",         database_path +"/frequencies_groups");
		System.setProperty("jing.chem.LJDatabase.pathName",                database_path +"/transport_groups");
		System.setProperty("jing.chem.TransportReferenceLibrary.pathName", database_path +"/transport_libraries");
		System.setProperty("jing.rxn.ReactionTemplateLibrary.pathName",    database_path +"/kinetics_groups");
		System.setProperty("jing.rxn.ReactionLibrary.pathName",            database_path +"/kinetics_libraries");
        System.setProperty("jing.chem.SolventLibrary.pathName",            database_path +"/thermo_libraries");
    }

    /**
     * Create a folder in the current (working) directory, deleting the
     * existing folder and its contents if desired.
     * @param name The name of the folder to create
     * @param deleteExisting true to delete the existing folder (and its contents!), false to preserve it
     */
    public static void createFolder(String name, boolean deleteExisting) {
        File folder = new File(name);
		if (deleteExisting)
            ChemParser.deleteDir(folder);
        if (!folder.exists())
            folder.mkdir();
    }

    public static void extractAndSetDatabasePath(String line) {
        StringTokenizer st = new StringTokenizer(line);
        String next = st.nextToken();
        String database_name = st.nextToken().trim();
        String workingDir = System.getProperty("RMG.workingDirectory");
        String database_path = workingDir + "/databases/" + database_name;
        setDatabasePaths(database_path);
    }

    /**
     * Echo the contents of the input condition file to the log.
     * @param inputfile The path of the input file
     */
    public static void logInputFile(String inputfile) throws FileNotFoundException, IOException {
        Logger.verbose("----------------------------------------------------------------------");
        Logger.verbose(" User input:");
        Logger.verbose("----------------------------------------------------------------------");
        Logger.verbose("");
        BufferedReader minireader = new BufferedReader(new FileReader(inputfile));
        String inputLine = minireader.readLine();
        while (inputLine != null) {
            Logger.verbose(inputLine);
            inputLine = minireader.readLine();
        }
        minireader.close();
        Logger.verbose("");
        Logger.verbose("----------------------------------------------------------------------");
        Logger.verbose("");
    }

}
/*********************************************************************
	File Path	: RMG\RMG\MainRMG.java
*********************************************************************/

