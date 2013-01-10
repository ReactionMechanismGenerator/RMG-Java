// //////////////////////////////////////////////////////////////////////////////
//
// RMG - Reaction Mechanism Generator
//
// Copyright (c) 2002-2011 Prof. William H. Green (whgreen@mit.edu) and the
// RMG Team (rmg_dev@mit.edu)
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// //////////////////////////////////////////////////////////////////////////////
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

// ----------------------------------------------------------------------------
// MainRMG.java
// ----------------------------------------------------------------------------
public class RMG {
    // ## configuration RMG::RMG
    public static void main(String[] args) {
        // Initialize the logger
        String logFileName = new File(System.getenv("RMG_JOB_OUTPUT"),
                "RMG.log").getPath();
        Logger.initialize(logFileName);
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
            createFolder(System.getProperty("RMG.jobOutputDir"), false);
            createFolder(System.getProperty("RMG.ChemkinOutputDir"), true);
            createFolder(System.getProperty("RMG.RestartDir"), false); // don't delete
            createFolder(System.getProperty("RMG.GATPFitDir"), true);
            createFolder(System.getProperty("RMG.ODESolverDir"), true);
            createFolder(System.getProperty("RMG.fameOutputDir"), true);
            createFolder(System.getProperty("RMG.frankieOutputDir"), true);
            createFolder(System.getProperty("RMG.jobScratchDir"), false);
            createFolder(System.getProperty("RMG.InChI_running_directory"),
                    true);
            createFolder(System.getProperty("RMG.PruningDir"), true);
            createFolder(System.getProperty("RMG.2DmolfilesDir"), true);
            createFolder(System.getProperty("RMG.3DmolfilesDir"), true);
            createFolder(System.getProperty("RMG.qmCalculationsDir"), false); // Preserving QM files between runs will
// speed things up considerably
            createFolder(System.getProperty("RMG.qmLibraryDir"), false); // don't delete
            // The only parameter should be the path to the condition file
            String inputfile = args[0];
            System.setProperty(
                    "jing.rxnSys.ReactionModelGenerator.conditionFile",
                    inputfile);
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
            Logger.info("Writing Final_Model.txt");
            writeFinalModel(rmg);
            // Delete remaining QM files if they were only meant tobe temporary
            if (!QMTP.keepQMfiles) {//keceli 20/Dec/12 There is a potential bug here 
                File qmFolder = new File(
                        System.getProperty("RMG.qmCalculationsDir"));
                if (qmFolder.exists()) {
                    ChemParser.deleteDir(qmFolder);
                }
            }
        } catch (Exception e) {
            // Any unhandled exception will land here
            // We assume that these were critical errors
            Logger.logStackTrace(e);
            Logger.critical(e.getMessage());
        }
        // Finish the logger
        Logger.finish();
        System.exit(0);
    }

    public static void writeFinalModel(ReactionModelGenerator rmg)
            throws IOException {
        CoreEdgeReactionModel cerm = (CoreEdgeReactionModel) rmg
                .getReactionModel();
        long tAtInitialization = Global.tAtInitialization;
        BufferedWriter file = new BufferedWriter(new FileWriter(
                "Final_Model.txt"));
        LinkedList speList = new LinkedList(rmg.getReactionModel()
                .getSpeciesSet());
        // 11/2/07 gmagoon: changing to loop over all reaction systems
        LinkedList rsList = rmg.getReactionSystemList();
        for (int i = 0; i < rsList.size(); i++) {
            ReactionSystem rs = (ReactionSystem) rsList.get(i);
            file.write("Reaction system " + (i + 1) + " of " + rsList.size()
                    + "\n");// 11/4/07 gmagoon: fixing "off by one" error
            file.write("\\\\\\\\\\\\\\\\\\\\\\\\\\\\    Concentration Profile Output    \\\\\\\\\\\\\\\\\\\\\\\\\\");
            file.write("\n" + rs.returnConcentrationProfile(speList) + "\n");
            if (rmg.getError()) {// svp
                file.write("\n");
                file.write("Upper Bounds:" + "\n");
                file.write(rs.printUpperBoundConcentrations(speList) + "\n");
                file.write("Lower Bounds:" + "\n");
                file.write(rs.printLowerBoundConcentrations(speList) + "\n");
            }
            file.write("\\\\\\\\\\\\\\\\\\\\\\\\\\\\    Mole Fraction Profile Output    \\\\\\\\\\\\\\\\\\\\\\\\\\");
            file.write("\n" + rs.returnMoleFractionProfile(speList) + "\n");
            file.write(rs.printOrderedReactions() + "\n");
            if (rmg.getSensitivity()) {// svp
                LinkedList importantSpecies = rmg.getSpeciesList();
                file.write("Sensitivity Analysis:");
                file.write(rs.printSensitivityCoefficients(speList,
                        importantSpecies) + "\n");
                file.write("\n");
                file.write(rs.printSensitivityToThermo(speList,
                        importantSpecies) + "\n");
                file.write("\n");
                file.write(rs.printMostUncertainReactions(speList,
                        importantSpecies) + "\n\n");
            }
            file.write(rs.returnReactionFlux() + "\n");
        }
        // Write how long the simulation took, and then close the file
        long end = System.currentTimeMillis();
        double min = (end - tAtInitialization) / 1.0E3 / 60.0;
        file.write(String.format("Running time: %.3f min", +min));
        file.close();
    }

    public static void globalInitializeSystemProperties() {
        // Set the default locale to US (since RMG is developed there)
        Locale.setDefault(Locale.US);
        // Get the value of the RMG environment variable; fail if not set!
        String workingDir = System.getenv("RMG");
        if (workingDir == null)
            throw new RuntimeException(
                    "The RMG environment variable is not defined.");
        System.setProperty("RMG.workingDirectory", workingDir);
        Logger.info(" Environment Variables:");
        Logger.info("RMG = " + workingDir);
        // Set the databases directory
        String databaseDir = System.getenv("RMG_DATABASES");
        if (databaseDir == null) {
            databaseDir = workingDir + "/databases";
            Logger.info("RMG_DATABASES = $RMG/databases (default)");
        } else {
            Logger.info(String.format(
                    "RMG_DATABASES = %s (environment variable)", databaseDir));
        }
        System.setProperty("RMG.databasesDirectory", databaseDir);
        // Set the default database. Will be over-ridden if a "Database:" line is found in the condition file.
        extractAndSetDatabasePath("Database: RMG_database");
        String source = "";
        // Set the job scratch dir
        String jobScratchDir = System.getenv("RMG_JOB_SCRATCH");
        if (jobScratchDir == null) {
            jobScratchDir = "."; // default of "scratch" may be nicer
            source = "default";
        } else
            source = "environment variable";
        Logger.info(String.format("RMG_JOB_SCRATCH = %s (%s)", jobScratchDir,
                source));
        System.setProperty("RMG.jobScratchDir", jobScratchDir);
        // Set the job output dir
        String jobOutputDir = System.getenv("RMG_JOB_OUTPUT");
        if (jobOutputDir == null) {
            jobOutputDir = ".";
            source = "default";
        } else
            source = "environment variable";
        Logger.info(String.format("RMG_JOB_OUTPUT = %s (%s)", jobOutputDir,
                source));
        System.setProperty("RMG.jobOutputDir", jobOutputDir);
        // Set the QM library directory
        String qmLibraryDir = System.getenv("RMG_QM_LIBRARY");
        if (qmLibraryDir == null) {
            qmLibraryDir = new File(jobOutputDir, "QMThermoLibrary").getPath();
            source = "default - relative to RMG_JOB_OUTPUT";
        } else
            source = "environment variable";
        Logger.info(String.format("RMG_QM_LIBRARY = %s (%s)", qmLibraryDir,
                source));
        System.setProperty("RMG.qmLibraryDir", qmLibraryDir);
        // Set the QM calculations directory
        String qmCalculationsDir = System.getenv("RMG_QM_CALCS");
        if (qmCalculationsDir == null) {
            qmCalculationsDir = new File(jobOutputDir, "QMfiles").getPath();
            source = "default - relative to RMG_JOB_OUTPUT";
        } else
            source = "environment variable";
        Logger.info(String.format("RMG_QM_CALCS = %s (%s)", qmCalculationsDir,
                source));
        System.setProperty("RMG.qmCalculationsDir", qmCalculationsDir);
        Logger.verbose(" Derived paths:");
        // Set the directory to save problematic Fame input/output in
        String fameOutputDir = new File(jobOutputDir, "fame").getPath();
        Logger.verbose("Fame errors directory = " + fameOutputDir);
        System.setProperty("RMG.fameOutputDir", fameOutputDir);
        // Set the directory to save problematic Frankie input/output in
        String frankieOutputDir = new File(jobOutputDir, "frankie").getPath();
        Logger.verbose("Frankie errors directory = " + frankieOutputDir);
        System.setProperty("RMG.frankieOutputDir", frankieOutputDir);
        // Set the directory to save problematic GATPFit input/output in
        String GATPFitDir = new File(jobOutputDir, "GATPFit").getPath();
        Logger.verbose("GATPFit errors directory = " + GATPFitDir);
        System.setProperty("RMG.GATPFitDir", GATPFitDir);
        // Set the directory to run the inchi executable in.
        System.setProperty("RMG.InChI_running_directory",
                new File(System.getProperty("RMG.jobScratchDir"), "InChI")
                        .getPath());
        Logger.verbose("InChI running directory = "
                + System.getProperty("RMG.InChI_running_directory"));
        // Set the directory to run the ODE solver in.
        System.setProperty("RMG.ODESolverDir",
                new File(System.getProperty("RMG.jobScratchDir"), "ODESolver")
                        .getPath());
        Logger.verbose("ODE Solver running directory = "
                + System.getProperty("RMG.ODESolverDir"));
        // 2D and 3D mol files for the RDKit portion of QM Thermo
        System.setProperty("RMG.2DmolfilesDir",
                new File(System.getProperty("RMG.jobScratchDir"), "2Dmolfiles")
                        .getPath());
        System.setProperty("RMG.3DmolfilesDir",
                new File(System.getProperty("RMG.jobScratchDir"), "3Dmolfiles")
                        .getPath());
        // output files
        // Set the directory to save the chemkin files in.
        System.setProperty("RMG.ChemkinOutputDir",
                new File(System.getProperty("RMG.jobOutputDir"), "chemkin")
                        .getPath());
        // Set the directory to save the Pruning files in.
        System.setProperty("RMG.PruningDir",
                new File(System.getProperty("RMG.jobOutputDir"), "Pruning")
                        .getPath());
        // Set the directory to save the Restart files in.
        System.setProperty("RMG.RestartDir",
                new File(System.getProperty("RMG.jobOutputDir"), "Restart")
                        .getPath());
    }

    public static void setDatabasePaths(String database_path) {
        // String database_path = workingDir + "/databases/" + name
        System.setProperty("jing.chem.ChemGraph.forbiddenStructureFile",
                database_path + "/ForbiddenStructures.txt");
        System.setProperty("jing.chem.ThermoGAGroupLibrary.pathName",
                database_path + "/thermo_groups");
        System.setProperty("jing.chem.ThermoReferenceLibrary.pathName",
                database_path + "/thermo_libraries");
        System.setProperty("jing.chem.FrequencyDatabase.pathName",
                database_path + "/frequencies_groups");
        System.setProperty("jing.chem.LJDatabase.pathName", database_path
                + "/transport_groups");
        System.setProperty("jing.chem.TransportReferenceLibrary.pathName",
                database_path + "/transport_libraries");
        System.setProperty("jing.rxn.ReactionTemplateLibrary.pathName",
                database_path + "/kinetics_groups");
        System.setProperty("jing.rxn.ReactionLibrary.pathName", database_path
                + "/kinetics_libraries");
        System.setProperty("jing.chem.SolventLibrary.pathName", database_path
                + "/thermo_libraries");
    }

    /**
     * Create a folder in the current (working) directory, deleting the existing folder and its contents if desired.
     * 
     * @param name
     *            The name of the folder to create
     * @param deleteExisting
     *            true to delete the existing folder (and its contents!), false to preserve it
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
        String databases_path = System.getProperty("RMG.databasesDirectory");
        String database_path = databases_path + "/" + database_name;
        setDatabasePaths(database_path);
    }

    /**
     * Echo the contents of the input condition file to the log.
     * 
     * @param inputfile
     *            The path of the input file
     */
    public static void logInputFile(String inputfile)
            throws FileNotFoundException, IOException {
        Logger.info("----------------------------------------------------------------------");
        Logger.info(" User input:");
        Logger.info("----------------------------------------------------------------------");
        Logger.info("");
        BufferedReader minireader = new BufferedReader(
                new FileReader(inputfile));
        String inputLine = minireader.readLine();
        while (inputLine != null) {
            Logger.info(inputLine);
            inputLine = minireader.readLine();
        }
        minireader.close();
        Logger.info("");
        Logger.info("----------------------------------------------------------------------");
        Logger.info("");
    }
}
/*********************************************************************
 * File Path : RMG\RMG\MainRMG.java
 *********************************************************************/
