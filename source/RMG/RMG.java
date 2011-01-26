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
            long tAtInitialization = System.currentTimeMillis();
            Global.tAtInitialization = tAtInitialization;
		
		    initializeSystemProperties(args[0]);
            ReactionModelGenerator rmg = new ReactionModelGenerator();

            // Generate the model!
            rmg.modelGeneration();

            CoreEdgeReactionModel cerm = (CoreEdgeReactionModel)rmg.getReactionModel();

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
                double min = (end - tAtInitialization) / 1E3 / 60;
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
       catch (Exception e) {
           // Any unhandled exception will land here
           // We assume that these were critial errors
           Logger.critical(e.getMessage());
       }

       // Finish the logger
       Logger.finish();
       
    }
    
    public static void globalInitializeSystemProperties() {
    	Locale.setDefault(Locale.US);
		String workingDir = System.getenv("RMG");
		 if (workingDir == null) {
			 System.err.println("Undefined environment variable RMG.");
			 System.exit(0);
		 }
		System.setProperty("RMG.workingDirectory", workingDir);
		String database_name= "RMG_database";
		String database_path = workingDir + "/databases/" + database_name;
		setDatabasePaths(database_path);
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

	public static void initializeSystemProperties(String inputfile) {
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

	     System.setProperty("jing.rxnSys.ReactionModelGenerator.conditionFile",inputfile);
		 try {
             String initialConditionFile = System.getProperty("jing.rxnSys.ReactionModelGenerator.conditionFile");
             if (initialConditionFile == null) {
                     System.out.println("undefined system property: jing.rxnSys.ReactionModelGenerator.conditionFile");
                     System.exit(0);
             }

			 SimpleDateFormat formatter = new SimpleDateFormat("yyyy-MM-dd'T'HH:mm:ss");
			 String formattedDate = formatter.format(new java.util.Date());
			 // to add timezone information, add 'Z' to the end of the SimpleDateFormat above and uncomment the following line
			 //formattedDate = formattedDate.substring(0, formattedDate.length()-2)  + ":" + formattedDate.substring(formattedDate.length()-2); // add the colon to the timezone
			 System.out.println("Current local time is: "+ formattedDate );
			 System.out.println("----------------------------------------------------------------------");
        	 System.out.println(" User input:");
             System.out.println("----------------------------------------------------------------------");

             //GJB: Added mini-reader to reprint input file in job output.
             FileReader tmpin = new FileReader(initialConditionFile);
             BufferedReader minireader = new BufferedReader(tmpin);
             try {
                 String inputLine = minireader.readLine();
                 while ( inputLine != null ) {
                	 System.out.println(inputLine);
                	 inputLine = minireader.readLine();
                 }
              }
             catch ( IOException error ) {
                 System.err.println( "Error reading condition file: " + error );
             }
             System.out.println("\n----------------------------------------------------------------------\n");
             minireader.close();
             
             FileReader in = new FileReader(initialConditionFile);
             BufferedReader reader = new BufferedReader(in);
             String line = ChemParser.readMeaningfulLine(reader, true);
			 //line = ChemParser.readMeaningfulLine(reader);

             if (line.startsWith("Database")){
                 extractAndSetDatabasePath(line);
			 }
	}
	catch (IOException e) {
		System.err.println("Error in read in reaction system initialization file!");
	}
  }

        public static void extractAndSetDatabasePath(String line) {
            StringTokenizer st = new StringTokenizer(line);
            String next = st.nextToken();
            String database_name = st.nextToken().trim();
            String workingDir = System.getProperty("RMG.workingDirectory");
            String database_path = workingDir + "/databases/" + database_name;
            setDatabasePaths(database_path);
        }

	
	/* this function is not used, has not been maintained, and is now out of date.
	private static void writeThermoFile() {
		// TODO Auto-generated method stub
		String thermoFile = "300.000  1000.000  5000.000 \n";
		thermoFile += "! neon added by pey (20/6/04) - used thermo for Ar\n";
		thermoFile += "Ne                120186Ne  1               G  0300.00   5000.00  1000.00      1\n";
		thermoFile += " 0.02500000E+02 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2\n";
		thermoFile += "-0.07453750E+04 0.04366001E+02 0.02500000E+02 0.00000000E+00 0.00000000E+00    3\n";
		thermoFile += " 0.00000000E+00 0.00000000E+00-0.07453750E+04 0.04366001E+02                   4\n";
		thermoFile += "N2                121286N   2               G  0300.00   5000.00  1000.00      1\n";
		thermoFile += " 0.02926640e+02 0.01487977e-01-0.05684761e-05 0.01009704e-08-0.06753351e-13    2\n";
		thermoFile += "-0.09227977e+04 0.05980528e+02 0.03298677e+02 0.01408240e-01-0.03963222e-04    3\n";
		thermoFile += " 0.05641515e-07-0.02444855e-10-0.01020900e+05 0.03950372e+02                   4\n";
		thermoFile += "Ar                120186Ar  1               G  0300.00   5000.00  1000.00      1\n";
		thermoFile += " 0.02500000e+02 0.00000000e+00 0.00000000e+00 0.00000000e+00 0.00000000e+00    2\n";
		thermoFile += "-0.07453750e+04 0.04366001e+02 0.02500000e+02 0.00000000e+00 0.00000000e+00    3\n";
		thermoFile += " 0.00000000e+00 0.00000000e+00-0.07453750e+04 0.04366001e+02                   4\n";
		thermoFile += "end\n";
		try {
			FileWriter fw = new FileWriter("chemkin/therm.dat");
			fw.write(thermoFile);
			fw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return;
	};
    */
}
/*********************************************************************
	File Path	: RMG\RMG\MainRMG.java
*********************************************************************/

